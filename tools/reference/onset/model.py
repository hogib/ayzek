"""The onset detector: a causal transformer over the near past, cross-attending
to a summary of the station's background.

    x (B, T, 4) ─ causal conv stem ─→ tokens (B, T/10, d)   one per 0.1 s
               ─ n_layers × [ sliding-window causal self-attention (ALiBi)
                              cross-attention to the station context
                              feed-forward ]
               ─ per token:  logit  "a P onset has arrived at or before this token"
                             dt     seconds since that onset

    ctx (B, Tc, 4) ─ same stem ─ K learned queries cross-attend ─→ (B, K, d)

Every output depends only on samples at or before the end of its token, so the
network can be run one token at a time (`stream.StreamingDetector`) with
exactly the outputs of a full forward pass. The tests hold it to that.

See docs/DESIGN.md for why each piece is here.
"""
from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F

from .config import DT_BIN_EDGES_S, FRESH_S, ModelConfig


# ---------------------------------------------------------------------------
# Stem
# ---------------------------------------------------------------------------

class CausalConv(nn.Module):
    """Conv1d padded on the left by `kernel - stride`.

    With that padding, output `i` covers input samples
    `[i*s - (k - s), (i + 1)*s - 1]`: it ends exactly on the last sample of its
    own stride. The more usual `k - 1` padding would leave each output blind to
    the newest `s - 1` samples, which is latency for nothing.
    """

    def __init__(self, cin: int, cout: int, kernel: int, stride: int):
        super().__init__()
        if kernel < stride:
            raise ValueError(f"kernel {kernel} < stride {stride} would skip samples")
        self.pad = kernel - stride
        self.conv = nn.Conv1d(cin, cout, kernel, stride=stride)

    def forward(self, x):
        return self.conv(F.pad(x, (self.pad, 0)))


class ChannelNorm(nn.Module):
    """LayerNorm over channels at each time step. Unlike BatchNorm or
    GroupNorm it never mixes time steps, so it cannot leak the future."""

    def __init__(self, channels: int):
        super().__init__()
        self.norm = nn.LayerNorm(channels)

    def forward(self, x):                       # (B, C, T)
        return self.norm(x.transpose(1, 2)).transpose(1, 2)


class Stem(nn.Module):
    """(B, T, C_in) samples -> (B, T / stride, d_model) tokens."""

    def __init__(self, cfg: ModelConfig):
        super().__init__()
        layers, cin = [], cfg.in_channels
        for k, s, c in cfg.stem:
            layers += [CausalConv(cin, c, k, s), ChannelNorm(c), nn.GELU()]
            cin = c
        self.net = nn.Sequential(*layers)
        self.proj = nn.Linear(cin, cfg.d_model)

    def forward(self, x):
        return self.proj(self.net(x.transpose(1, 2)).transpose(1, 2))


# ---------------------------------------------------------------------------
# Attention
# ---------------------------------------------------------------------------

def alibi_slopes(n_heads: int) -> torch.Tensor:
    """Per-head distance penalties (Press et al. 2022). With 4 heads they run
    from 1/4 to 1/256 per token, so one head is nearly local and one sees the
    whole window evenly."""
    return torch.tensor([2.0 ** (-8.0 * (i + 1) / n_heads) for i in range(n_heads)])


def sliding_alibi_bias(n_tokens: int, window: int, slopes: torch.Tensor) -> torch.Tensor:
    """(H, T, T) additive bias: -slope * distance inside the causal window,
    -inf outside it. The diagonal is always inside, so no row is empty."""
    i = torch.arange(n_tokens, device=slopes.device)
    dist = (i[:, None] - i[None, :]).float()
    bias = -slopes[:, None, None] * dist
    outside = (dist < 0) | (dist >= window)
    return bias.masked_fill(outside, float("-inf"))


class Attention(nn.Module):
    """Multi-head attention with separate query and key/value inputs, and the
    projections exposed so the streaming path can cache keys and values."""

    def __init__(self, d: int, n_heads: int, dropout: float):
        super().__init__()
        self.h, self.dk = n_heads, d // n_heads
        self.q = nn.Linear(d, d)
        self.kv = nn.Linear(d, 2 * d)
        self.out = nn.Linear(d, d)
        self.dropout = dropout

    def split(self, x):                          # (B, T, d) -> (B, H, T, dk)
        B, T, _ = x.shape
        return x.view(B, T, self.h, self.dk).transpose(1, 2)

    def project_kv(self, mem):
        k, v = self.kv(mem).chunk(2, dim=-1)
        return self.split(k), self.split(v)

    def attend(self, q, k, v, bias=None):
        if bias is not None:
            bias = bias.to(q.dtype)
        o = F.scaled_dot_product_attention(
            q, k, v, attn_mask=bias,
            dropout_p=self.dropout if self.training else 0.0)
        B, H, T, dk = o.shape
        return self.out(o.transpose(1, 2).reshape(B, T, H * dk))

    def forward(self, x, mem, bias=None):
        k, v = self.project_kv(mem)
        return self.attend(self.split(self.q(x)), k, v, bias)


class FeedForward(nn.Sequential):
    def __init__(self, d: int, mult: int, dropout: float):
        super().__init__(nn.Linear(d, d * mult), nn.GELU(), nn.Dropout(dropout),
                         nn.Linear(d * mult, d))


class Block(nn.Module):
    """Pre-LN: causal self-attention, cross-attention to context, feed-forward."""

    def __init__(self, cfg: ModelConfig):
        super().__init__()
        d = cfg.d_model
        self.ln_sa, self.ln_ca, self.ln_ff = nn.LayerNorm(d), nn.LayerNorm(d), nn.LayerNorm(d)
        self.sa = Attention(d, cfg.n_heads, cfg.dropout)
        self.ca = Attention(d, cfg.n_heads, cfg.dropout)
        self.ff = FeedForward(d, cfg.ffn_mult, cfg.dropout)
        self.drop = nn.Dropout(cfg.dropout)

    def forward(self, x, ctx, bias):
        h = self.ln_sa(x)
        x = x + self.drop(self.sa(h, h, bias))
        x = x + self.drop(self.ca(self.ln_ca(x), ctx))
        return x + self.drop(self.ff(self.ln_ff(x)))


# ---------------------------------------------------------------------------
# Station context
# ---------------------------------------------------------------------------

class ContextEncoder(nn.Module):
    """A minute or two of the station's quiet signal -> K summary tokens.

    K learned queries cross-attend to the stem tokens of the context
    (Perceiver-style), then attend among themselves. There is no positional
    term: the context is a description of the background, not a sequence to
    be read in order. Tokens that span a gap are masked out.
    """

    def __init__(self, cfg: ModelConfig, stem: Stem):
        super().__init__()
        d = cfg.d_model
        self.stem = stem
        self.queries = nn.Parameter(torch.randn(cfg.context_tokens, d) * 0.02)
        self.ln_mem, self.ln_q = nn.LayerNorm(d), nn.LayerNorm(d)
        self.cross = Attention(d, cfg.n_heads, cfg.dropout)
        self.ln_ff = nn.LayerNorm(d)
        self.ff = FeedForward(d, cfg.ffn_mult, cfg.dropout)
        self.self_layers = nn.ModuleList(
            nn.ModuleDict({"ln_sa": nn.LayerNorm(d),
                           "sa": Attention(d, cfg.n_heads, cfg.dropout),
                           "ln_ff": nn.LayerNorm(d),
                           "ff": FeedForward(d, cfg.ffn_mult, cfg.dropout)})
            for _ in range(cfg.context_layers))
        self.ln_out = nn.LayerNorm(d)
        self.stride = cfg.stride

    def forward(self, ctx):                     # (B, Tc, 4) -> (B, K, d)
        mem = self.ln_mem(self.stem(ctx))
        n = mem.shape[1]
        gap = ctx[:, : n * self.stride, 3].reshape(ctx.shape[0], n, self.stride).amax(-1) > 0
        # A fully gapped context would leave softmax nothing to attend to;
        # callers mark those examples has_ctx=False, and unmasking here keeps
        # the arithmetic finite for them.
        gap = gap & ~gap.all(dim=1, keepdim=True)
        bias = torch.zeros(gap.shape, device=ctx.device).masked_fill(gap, float("-inf"))
        bias = bias[:, None, None, :]           # (B, 1, 1, Tc_tokens)

        q = self.queries.expand(ctx.shape[0], -1, -1)
        z = q + self.cross(self.ln_q(q), mem, bias)
        z = z + self.ff(self.ln_ff(z))
        for layer in self.self_layers:
            h = layer["ln_sa"](z)
            z = z + layer["sa"](h, h)
            z = z + layer["ff"](layer["ln_ff"](z))
        return self.ln_out(z)


# ---------------------------------------------------------------------------
# The detector
# ---------------------------------------------------------------------------

class OnsetDetector(nn.Module):
    def __init__(self, cfg: ModelConfig):
        super().__init__()
        self.cfg = cfg
        d = cfg.d_model
        self.stem = Stem(cfg)
        self.context = ContextEncoder(cfg, self.stem)
        # Used in place of an encoded context when there is none: a station
        # that has just come up, a gapped background, or context dropout.
        self.null_context = nn.Parameter(torch.randn(cfg.context_tokens, d) * 0.02)
        self.blocks = nn.ModuleList(Block(cfg) for _ in range(cfg.n_layers))
        self.ln_out = nn.LayerNorm(d)
        self.head = nn.Linear(d, 2)
        if cfg.geometry:
            # log distance and its log-variance.
            self.geo_head = nn.Linear(d, 2)
        if cfg.dt_bins:
            self.dt_head = nn.Linear(d, len(DT_BIN_EDGES_S) + 1)
        self.register_buffer("slopes", alibi_slopes(cfg.n_heads), persistent=False)

    def encode_context(self, ctx, has_ctx, batch_size: int):
        """(B, K, d) context tokens; the null context where `has_ctx` is False."""
        null = self.null_context.expand(batch_size, -1, -1)
        if ctx is None or has_ctx is None or not bool(has_ctx.any()):
            return null
        enc = self.context(ctx)
        return torch.where(has_ctx[:, None, None], enc, null)

    def readout(self, h):
        """Final hidden state -> (logit, dt seconds).

        Always in fp32. A trained detector's logits sit around 5–12, where bf16
        resolves only ~0.03–0.06, so under autocast the probabilities near 1
        (exactly where the operating thresholds are) come out quantized.
        """
        with torch.autocast(h.device.type, enabled=False):
            out = self.head(self.ln_out(h.float()))
        return out[..., 0], self.cfg.max_dt_s * torch.sigmoid(out[..., 1])

    def dt_bins(self, h):
        """Final hidden state -> logits over the dt bins (DT_BIN_EDGES_S), in
        fp32 like `readout`. The dt regression answers "how long ago"; these
        say how sure the model is that the onset is fresh, which a single
        number cannot: an uncertain regression settles between the values,
        where a level trigger reads it as no."""
        with torch.autocast(h.device.type, enabled=False):
            return self.dt_head(self.ln_out(h.float()))

    def geometry(self, h):
        """Final hidden state -> how far away the event is, as seen from this
        station: a Gaussian in log km.

        Only meaningful after P. Before S the model can only infer distance from
        the P wave itself; once S is inside its lookback it can in effect read
        the S-P time, and the stated uncertainty should shrink to match.
        """
        with torch.autocast(h.device.type, enabled=False):
            g = self.geo_head(self.ln_out(h.float()))
        min_lv = 2.0 * math.log(self.cfg.geo_min_sd)
        return {"log_dist": g[..., 0], "log_dist_var": g[..., 1].clamp(min_lv, 6.0)}

    def forward(self, x, ctx=None, has_ctx=None):
        """
        Args:
            x: (B, T, 4) conditioned samples, T a multiple of `cfg.stride`.
            ctx: (B, Tc, 4) conditioned station context, or None.
            has_ctx: (B,) bool, which rows of `ctx` are real.

        Returns:
            dict with `logit` and `dt`, each (B, T / stride).
        """
        if x.shape[1] % self.cfg.stride:
            raise ValueError(f"sequence length {x.shape[1]} is not a multiple of "
                             f"the token stride {self.cfg.stride}")
        h = self.stem(x)
        c = self.encode_context(ctx, has_ctx, x.shape[0])
        bias = sliding_alibi_bias(h.shape[1], self.cfg.window_tokens, self.slopes)
        for block in self.blocks:
            h = block(h, c, bias[None])
        logit, dt = self.readout(h)
        out = {"logit": logit, "dt": dt}
        if self.cfg.geometry:
            out.update(self.geometry(h))
        if self.cfg.dt_bins:
            out["dt_logits"] = self.dt_bins(h)
        return out


def dt_bin_targets(dt: torch.Tensor) -> torch.Tensor:
    """dt targets in seconds -> bin indices on DT_BIN_EDGES_S."""
    return torch.bucketize(dt, torch.tensor(DT_BIN_EDGES_S, device=dt.device), right=True)


def fresh_probability(dt_logits: torch.Tensor) -> torch.Tensor:
    """P(dt < FRESH_S): the probability mass of the bins below FRESH_S."""
    k = sum(e <= FRESH_S for e in DT_BIN_EDGES_S)
    return torch.softmax(dt_logits.float(), dim=-1)[..., :k].sum(-1)


def count_parameters(model: nn.Module) -> int:
    return sum(p.numel() for p in model.parameters())
