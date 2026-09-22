"""The P/S picker, `wave` arm: strided convolutions into a BiLSTM, dense output.

Copied from sphase (`src/sphase/model.py`, at fd45c16). Only the `wave` arm is
kept, which is the one ayzek ships; sphase's `ram`, `beta` and dual arms are
left out, and with them its RAM-image module. Parameter names are unchanged, so
sphase checkpoints of the `wave` arm load strictly.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F


class WaveBranch(nn.Module):
    """Strided convolutions down to chunk resolution, then a BiLSTM."""

    def __init__(self, in_ch: int = 3, width: int = 32, out_dim: int = 64,
                 n_chunks: int = 250):
        super().__init__()
        self.n_chunks = n_chunks
        self.stem = nn.Sequential(
            nn.Conv1d(in_ch, width, 9, stride=2, padding=4), nn.BatchNorm1d(width), nn.ReLU(),
            nn.Conv1d(width, width * 2, 9, stride=2, padding=4), nn.BatchNorm1d(width * 2), nn.ReLU(),
            nn.Conv1d(width * 2, width * 2, 7, stride=2, padding=3), nn.BatchNorm1d(width * 2), nn.ReLU(),
        )
        self.lstm = nn.LSTM(width * 2, out_dim // 2, batch_first=True,
                            bidirectional=True, num_layers=1)
        self.attn = nn.MultiheadAttention(out_dim, num_heads=4, batch_first=True)
        self.norm = nn.LayerNorm(out_dim)

    def forward(self, x: torch.Tensor) -> torch.Tensor:      # (B, 3, m) -> (B, n, D)
        h = self.stem(x)
        h = F.adaptive_avg_pool1d(h, self.n_chunks).transpose(1, 2)
        h, _ = self.lstm(h)
        a, _ = self.attn(h, h, h, need_weights=False)
        return self.norm(h + a)


class PhasePicker(nn.Module):
    """`sphase.model.PhasePicker(arm="wave")`: one distribution over {noise, P, S} per chunk."""

    def __init__(self, n_chunks: int = 250, dim: int = 64):
        super().__init__()
        self.arm, self.n_chunks = "wave", n_chunks
        self.wave = WaveBranch(out_dim=dim, n_chunks=n_chunks)
        self.head = nn.Sequential(nn.Linear(dim, dim), nn.ReLU(), nn.Linear(dim, 3))

    def forward(self, x: torch.Tensor) -> torch.Tensor:      # (B, 3, m) -> (B, n, 3) logits
        return self.head(self.wave(x))
