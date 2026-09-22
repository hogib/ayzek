"""The detector as trained: one waveform branch and a binary head.

Copied from archive_pipeline (`products/detector.py`, at db13814), which loads
the training checkpoints for scoring. Parameter names match the training model's
(`b1`, `p1`, `w1`, `w2`, `head`) because the checkpoints are keyed on them.
"""
import torch
import torch.nn as nn


class ConvSeqBranch(nn.Module):
    """Strided 1D convolutions, then BiLSTM and self-attention.

    The convolutions reduce the sequence about 8x before the recurrent layer,
    which also makes the attention that follows ~64x cheaper. `use_lstm=False`
    stops after the convolutions and mean-pools.
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
        """Encodes `(batch, time, in_dim)` into `(batch, out_dim)`."""
        h = self.conv(x.transpose(1, 2)).transpose(1, 2)
        if not self.use_lstm:
            return h.mean(dim=1)
        h, _ = self.lstm(h)
        a, _ = self.attn(h, h, h)
        return self.norm(h + a).mean(dim=1)


class LSTMAttentionBranch(nn.Module):
    """BiLSTM over raw samples, then self-attention. The pre-convolution design."""

    def __init__(self, in_dim, hidden=64, layers=1, heads=4, dropout=0.2):
        super().__init__()
        self.lstm = nn.LSTM(in_dim, hidden, num_layers=layers, batch_first=True,
                            bidirectional=True,
                            dropout=dropout if layers > 1 else 0.0)
        d = hidden * 2
        self.attn = nn.MultiheadAttention(d, heads, dropout=dropout,
                                          batch_first=True)
        self.norm = nn.LayerNorm(d)
        self.out_dim = d

    def forward(self, x):
        """Encodes `(batch, time, in_dim)` into `(batch, out_dim)`."""
        h, _ = self.lstm(x)
        a, _ = self.attn(h, h, h)
        return self.norm(h + a).mean(dim=1)


class WaveformDetector(nn.Module):
    """One waveform branch and a binary head: `--channels 1d`, `--fusion linear`.

    Parameter names match the training model's exactly -- `b1`, `p1`, `w1`,
    `w2`, `head` -- because the checkpoints are keyed on them.
    """

    def __init__(self, seq_dim=3, hidden=48, fusion_dim=96, dropout=0.4,
                 branch1d="cnn-lstm", lstm_layers=1, lstm_heads=4):
        """Builds the network.

        Args:
            seq_dim: Channels per timestep; 3 for Z/N/E.
            hidden: LSTM hidden size per direction.
            fusion_dim: Width the branch is projected to, and the head's width.
            dropout: Inactive at eval, but it sizes nothing, so it is carried
                only so a model built here can also be trained elsewhere.
            branch1d: `cnn-lstm`, `cnn` or `lstm`.
            lstm_layers: Stacked LSTM layers.
            lstm_heads: Attention heads; must divide `hidden * 2`.

        Raises:
            ValueError: On an unknown `branch1d`.
        """
        super().__init__()
        if branch1d == "lstm":
            self.b1 = LSTMAttentionBranch(seq_dim, hidden=hidden,
                                          layers=lstm_layers, heads=lstm_heads,
                                          dropout=dropout)
        elif branch1d in ("cnn", "cnn-lstm"):
            self.b1 = ConvSeqBranch(seq_dim, hidden=hidden, layers=lstm_layers,
                                    heads=lstm_heads, dropout=dropout,
                                    use_lstm=(branch1d == "cnn-lstm"))
        else:
            raise ValueError(
                f"branch1d must be 'lstm', 'cnn' or 'cnn-lstm', got {branch1d!r}")
        self.p1 = nn.Linear(self.b1.out_dim, fusion_dim)
        self.w1 = nn.Parameter(torch.tensor(1.0))
        self.w2 = nn.Parameter(torch.tensor(1.0))
        self.head = nn.Sequential(
            nn.LayerNorm(fusion_dim),
            nn.Dropout(dropout),
            nn.Linear(fusion_dim, fusion_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(fusion_dim, 1),
        )

    def forward(self, seq):
        """Raw logit per window, shape `(batch, 1)`."""
        return self.head(self.w1 * self.p1(self.b1(seq)))
