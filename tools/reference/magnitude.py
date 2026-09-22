"""The magnitude regressor's network: the dual-channel trunk and its branches.

Copied from cascade_impl (`blocks.py` and `models.py`, at 1e4f6bf), the
project that trains it. The shipped checkpoints load strictly into
`DualChannelNet(3, 3, aux_dim=0, hidden=64, fusion_dim=128, dropout=0.3,
channels="1d+2d", n_classes=1, squeeze_output=True, branch1d="lstm")`.
"""
import torch
import torch.nn as nn


class LSTMAttentionBranch(nn.Module):
    """LSTM for long-range order, then multi-head self-attention to weight steps."""

    def __init__(self, in_dim, hidden=64, layers=1, heads=4, dropout=0.2):
        """Initializes the LSTM+attention branch.

        Args:
            in_dim: Size of each step's input feature vector.
            hidden: LSTM hidden size per direction; the bidirectional output
                (and attention/LayerNorm width) is ``hidden * 2``.
            layers: Number of stacked LSTM layers.
            heads: Number of attention heads. Must divide ``hidden * 2``.
            dropout: Dropout used inside the LSTM (when ``layers > 1``) and
                the attention module.
        """
        super().__init__()
        self.lstm = nn.LSTM(in_dim, hidden, num_layers=layers, batch_first=True,
                            bidirectional=True,
                            dropout=dropout if layers > 1 else 0.0)
        d = hidden * 2
        self.attn = nn.MultiheadAttention(d, heads, dropout=dropout, batch_first=True)
        self.norm = nn.LayerNorm(d)
        self.out_dim = d

    def forward(self, x):
        """Encodes a sequence into one pooled embedding.

        Args:
            x: Input sequence, shape (batch, time, in_dim).

        Returns:
            Tensor of shape (batch, out_dim), the time-mean of the
            attention+residual output.
        """
        h, _ = self.lstm(x)
        a, _ = self.attn(h, h, h)
        h = self.norm(h + a)             # residual, as in the transformer block
        return h.mean(dim=1)             # pool over time


class ConvSeqBranch(nn.Module):
    """1D CNN over the raw waveform, optionally followed by BiLSTM + attention.

    `LSTMAttentionBranch` feeds 600 raw 100 Hz samples straight into an LSTM,
    so the only local structure available to it is whatever recurrence can
    accumulate one 10 ms sample at a time. The detectors this project compares
    itself against do the opposite: EQTransformer is CNN -> BiLSTM -> attention
    and PhaseNet is a U-Net of 1D convolutions, both extracting local waveform
    features convolutionally before any recurrence.

    This branch adds that missing front end. Strided convolutions reduce the
    sequence roughly 8x before the recurrent layer, which also makes the
    self-attention that follows ~64x cheaper (its cost is quadratic in
    sequence length).

    `use_lstm=False` stops after the convolutions and mean-pools, isolating
    whether the recurrence contributes anything once local features exist.

    Args:
        in_dim: Channels per timestep of the input sequence (3 for Z/N/E).
        hidden: LSTM hidden size per direction. Output width is ``hidden * 2``
            with an LSTM, or ``conv_width`` without one.
        layers: Stacked LSTM layers.
        heads: Attention heads. Must divide ``hidden * 2``.
        dropout: Dropout after the convolution stages and inside attention.
        conv_width: Channel width of the final convolution stage.
    """

    def __init__(self, in_dim, hidden=64, layers=1, heads=4, dropout=0.2,
                 use_lstm=True, conv_width=96):
        super().__init__()
        self.use_lstm = use_lstm

        def stage(cin, cout, k, s):
            return nn.Sequential(
                nn.Conv1d(cin, cout, kernel_size=k, stride=s, padding=k // 2,
                          bias=False),
                nn.BatchNorm1d(cout), nn.GELU())

        # 600 -> 300 -> 150 -> 75 at 100 Hz. Kernels stay wide enough at the
        # first stage (70 ms) to see an arrival's onset rather than one cycle.
        self.conv = nn.Sequential(
            stage(in_dim, conv_width // 4, 7, 2),
            stage(conv_width // 4, conv_width // 2, 5, 2),
            nn.Dropout(dropout),
            stage(conv_width // 2, conv_width, 5, 2),
        )

        if use_lstm:
            self.lstm = nn.LSTM(conv_width, hidden, num_layers=layers,
                                batch_first=True, bidirectional=True,
                                dropout=dropout if layers > 1 else 0.0)
            d = hidden * 2
            self.attn = nn.MultiheadAttention(d, heads, dropout=dropout,
                                              batch_first=True)
            self.norm = nn.LayerNorm(d)
            self.out_dim = d
        else:
            self.out_dim = conv_width

    def forward(self, x):
        """Encodes a raw waveform sequence into one pooled embedding.

        Args:
            x: Input sequence, shape (batch, time, in_dim).

        Returns:
            Tensor of shape (batch, out_dim), mean-pooled over time.
        """
        h = self.conv(x.transpose(1, 2)).transpose(1, 2)   # (B, T', conv_width)
        if not self.use_lstm:
            return h.mean(dim=1)
        h, _ = self.lstm(h)
        a, _ = self.attn(h, h, h)
        h = self.norm(h + a)
        return h.mean(dim=1)


class CNNBranch(nn.Module):
    """Compact CNN over the 2D channel -- the log-power spectrogram.

    The images are small (33x38 for the detector's n_fft 64, 129x10 for the
    regressor's n_fft 256), so a 4-stage ResNet would be heavily
    over-provisioned here.
    """

    def __init__(self, in_channels=3, width=32, dropout=0.2):
        """Initializes the compact image branch.

        Args:
            in_channels: Number of input image channels.
            width: Base channel width; the three conv stages use
                ``width``, ``width*2``, ``width*4`` channels.
            dropout: Dropout2d probability applied after the second stage.
        """
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(in_channels, width, 3, padding=1, bias=False),
            nn.BatchNorm2d(width), nn.GELU(),
            nn.Conv2d(width, width * 2, 3, stride=2, padding=1, bias=False),
            nn.BatchNorm2d(width * 2), nn.GELU(),
            nn.Dropout2d(dropout),
            nn.Conv2d(width * 2, width * 4, 3, stride=2, padding=1, bias=False),
            nn.BatchNorm2d(width * 4), nn.GELU(),
            nn.AdaptiveAvgPool2d(1),
        )
        self.out_dim = width * 4

    def forward(self, x):
        """Encodes an image into one pooled embedding.

        Args:
            x: Input image batch, shape (batch, in_channels, height, width).

        Returns:
            Tensor of shape (batch, out_dim).
        """
        return torch.flatten(self.net(x), 1)


class GatedFusion(nn.Module):
    """
    Per-example gate deciding how much to trust each branch, replacing a
    fixed pair of scalars (a*F1 + b*F2, same for every example) with
    g(x)*F1 + (1-g(x))*F2, where g = sigmoid(MLP([F1, F2])) is conditioned on
    both branches' own features for THIS example.
    """

    def __init__(self, dim, hidden=None, dropout=0.1):
        """Initializes the gated fusion module.

        Args:
            dim: Width of each of the two branch embeddings being fused
                (they must match).
            hidden: Hidden width of the gating MLP. Defaults to
                ``max(8, dim // 2)`` when None.
            dropout: Dropout probability inside the gating MLP.
        """
        super().__init__()
        hidden = hidden or max(8, dim // 2)
        self.net = nn.Sequential(
            nn.Linear(dim * 2, hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, 1),
        )

    def forward(self, f1, f2):
        """Fuses two branch embeddings with a learned per-example gate.

        Args:
            f1: First branch's embedding, shape (batch, dim).
            f2: Second branch's embedding, shape (batch, dim).

        Returns:
            Tuple of (fused, gate): ``fused`` has shape (batch, dim);
            ``gate`` has shape (batch, 1) and is in (0, 1), where values
            near 1 favor ``f1`` and values near 0 favor ``f2``.
        """
        g = torch.sigmoid(self.net(torch.cat([f1, f2], dim=1)))
        return g * f1 + (1.0 - g) * f2, g


class DualChannelTrunk(nn.Module):
    """Builds the 1D/2D/aux branches and fuses them into one vector of width
    `fused_dim`. Shared by every dual-channel model in this repo; subclasses
    attach whatever head(s) their task needs.

    `channels` ablates which branches are active: "all", "1d", "2d", "aux",
    "1d+aux", "2d+aux" (aux-only variants require `aux_dim > 0`), and "1d+2d"
    (both waveform branches with the aux vector withheld).
    `fusion="linear"` is a*F1+b*F2 with learned scalars (the paper's
    default); `fusion="gate"` is a per-example gate (`blocks.GatedFusion`),
    only meaningful when both 1D and 2D are active.
    """

    def __init__(self, seq_dim, img_channels, aux_dim=0, hidden=64, fusion_dim=128,
                dropout=0.3, channels="all", fusion="linear",
                lstm_layers=1, lstm_heads=4, branch1d="lstm"):
        """Initializes the 1D/2D/aux branches and fusion.

        Args:
            seq_dim: Per-step feature width of the 1D sequence input.
            img_channels: Number of channels of the 2D image input.
            aux_dim: Width of an auxiliary scalar vector concatenated after
                fusion. 0 disables the aux branch.
            hidden: LSTM hidden size (per direction) for the 1D branch.
            fusion_dim: Common width both branches are projected to before
                fusion, and the fused output's width.
            dropout: Dropout used throughout the branches and fusion.
            channels: Which branches are active -- "all", "1d", "2d", "aux",
                "1d+aux", "2d+aux" (aux-only variants require `aux_dim > 0`),
                or "1d+2d" (both waveform branches, no aux -- what a cascade
                can actually supply at run time).
            fusion: "linear" (a*F1+b*F2 with learned scalars) or "gate" (a
                per-example gate, `blocks.GatedFusion`); "gate" only
                takes effect when both 1D and 2D are active.
            lstm_layers: Number of stacked LSTM layers in the 1D branch.
            lstm_heads: Number of attention heads in the 1D branch.
            branch1d: Architecture of the 1D branch. "lstm" (default) is
                `LSTMAttentionBranch`, which reads raw samples directly and is
                what every existing result used. "cnn-lstm" prepends a strided
                1D convolutional encoder (`ConvSeqBranch`), matching the
                CNN-then-recurrence order EQTransformer and PhaseNet use;
                "cnn" keeps the convolutions and drops the recurrence.

        Raises:
            ValueError: If `channels` disables every branch, or `fusion` is
                not "linear" or "gate".
        """
        super().__init__()
        self.channels = channels
        self.fusion = fusion
        self.aux_dim = aux_dim
        self.use_1d = channels in ("all", "1d", "1d+aux", "1d+2d")
        self.use_2d = channels in ("all", "2d", "2d+aux", "1d+2d")
        # "1d+2d" is deliberately absent here: both waveform branches, no aux.
        # It exists because a deployable cascade cannot supply the aux vector.
        # `log_distance` is the epicentral distance to a CATALOGUED hypocentre,
        # and a window the detector just flagged has no catalogue entry, so the
        # one configuration an operational stage 2 needs -- everything the
        # waveform gives and nothing it does not -- had no name.
        self.use_aux = aux_dim > 0 and channels in ("all", "aux", "1d+aux", "2d+aux")
        if not (self.use_1d or self.use_2d or self.use_aux):
            raise ValueError(f"--channels {channels} disables every branch")
        if fusion not in ("linear", "gate"):
            raise ValueError(f"--fusion must be 'linear' or 'gate', got {fusion!r}")

        if branch1d not in ("lstm", "cnn", "cnn-lstm"):
            raise ValueError(
                f"branch1d must be 'lstm', 'cnn' or 'cnn-lstm', got {branch1d!r}")
        self.branch1d = branch1d

        if self.use_1d:
            if branch1d == "lstm":
                self.b1 = LSTMAttentionBranch(seq_dim, hidden=hidden,
                                              layers=lstm_layers,
                                              heads=lstm_heads, dropout=dropout)
            else:
                self.b1 = ConvSeqBranch(seq_dim, hidden=hidden, layers=lstm_layers,
                                        heads=lstm_heads, dropout=dropout,
                                        use_lstm=(branch1d == "cnn-lstm"))
            self.p1 = nn.Linear(self.b1.out_dim, fusion_dim)
        if self.use_2d:
            self.b2 = CNNBranch(img_channels, dropout=dropout)
            self.p2 = nn.Linear(self.b2.out_dim, fusion_dim)

        self.both = self.use_1d and self.use_2d
        if self.both and fusion == "gate":
            self.gated_fusion = GatedFusion(fusion_dim)
        else:
            # Learned fusion weights (a, b in the paper's notation). Also
            # used, harmlessly, as a global rescale in single-branch
            # ablations -- the optimizer settles it near 1 since there is
            # nothing to balance it against.
            self.w1 = nn.Parameter(torch.tensor(1.0))
            self.w2 = nn.Parameter(torch.tensor(1.0))

        self.fused_dim = (fusion_dim if (self.use_1d or self.use_2d) else 0) + \
                         (aux_dim if self.use_aux else 0)
        self.last_gate = None

    def _fuse(self, seq, img, aux):
        """Runs the active branches and fuses them into one vector.

        Args:
            seq: 1D sequence input, shape (batch, time, seq_dim). Unused if
                `use_1d` is False.
            img: 2D image input, shape (batch, img_channels, height, width).
                Unused if `use_2d` is False.
            aux: Auxiliary scalar input, shape (batch, aux_dim). Unused if
                `use_aux` is False.

        Returns:
            Tensor of shape (batch, fused_dim). Also sets `self.last_gate`
            to the per-example gate (batch, 1) when `fusion="gate"` and both
            branches are active, else None.
        """
        self.last_gate = None
        feats = []
        fused = None
        if self.both:
            f1 = self.p1(self.b1(seq))
            f2 = self.p2(self.b2(img))
            if self.fusion == "gate":
                fused, self.last_gate = self.gated_fusion(f1, f2)
            else:
                fused = self.w1 * f1 + self.w2 * f2
        elif self.use_1d:
            fused = self.w1 * self.p1(self.b1(seq))
        elif self.use_2d:
            fused = self.w2 * self.p2(self.b2(img))
        if fused is not None:
            feats.append(fused)
        if self.use_aux:
            feats.append(aux)
        return torch.cat(feats, dim=1)


class DualChannelNet(DualChannelTrunk):
    """`DualChannelTrunk` plus a single Sequential head of width `n_classes`
    (1 for binary/regression, 3+ for multiclass). `squeeze_output` matches
    each caller's existing convention for a single-logit head."""

    def __init__(self, seq_dim, img_channels, aux_dim=0, hidden=64, fusion_dim=128,
                dropout=0.3, channels="all", fusion="linear",
                lstm_layers=1, lstm_heads=4, n_classes=1, squeeze_output=False,
                branch1d="lstm"):
        """Initializes the trunk (see `DualChannelTrunk.__init__`) plus a head.

        Args:
            seq_dim: See `DualChannelTrunk.__init__`.
            img_channels: See `DualChannelTrunk.__init__`.
            aux_dim: See `DualChannelTrunk.__init__`.
            hidden: See `DualChannelTrunk.__init__`.
            fusion_dim: See `DualChannelTrunk.__init__`.
            dropout: See `DualChannelTrunk.__init__`.
            channels: See `DualChannelTrunk.__init__`.
            fusion: See `DualChannelTrunk.__init__`.
            lstm_layers: See `DualChannelTrunk.__init__`.
            lstm_heads: See `DualChannelTrunk.__init__`.
            n_classes: Width of the head's output layer -- 1 for
                binary/regression, 3+ for multiclass.
            squeeze_output: If True, squeezes the last dimension off the
                output (for a single-logit head returned as shape (batch,)
                rather than (batch, 1)).
        """
        super().__init__(seq_dim, img_channels, aux_dim=aux_dim, hidden=hidden,
                         fusion_dim=fusion_dim, dropout=dropout, channels=channels,
                         fusion=fusion, lstm_layers=lstm_layers, lstm_heads=lstm_heads,
                         branch1d=branch1d)
        self.squeeze_output = squeeze_output
        self.head = nn.Sequential(
            nn.LayerNorm(self.fused_dim),
            nn.Dropout(dropout),
            nn.Linear(self.fused_dim, fusion_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(fusion_dim, n_classes),
        )

    def forward(self, seq, img, aux=None):
        """Fuses the branches and applies the head.

        Args:
            seq: See `DualChannelTrunk._fuse`.
            img: See `DualChannelTrunk._fuse`.
            aux: See `DualChannelTrunk._fuse`.

        Returns:
            Tensor of shape (batch, n_classes), or (batch,) if
            `squeeze_output` and `n_classes == 1`.
        """
        out = self.head(self._fuse(seq, img, aux))
        return out.squeeze(-1) if self.squeeze_output else out
