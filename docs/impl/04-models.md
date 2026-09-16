# 04 · Models

`src/nn.cpp` holds the layers and `src/models.cpp` composes them into the two
trained networks. Each layer matches its PyTorch module in eval mode and reads
the same `state_dict` names.

## Layers

| layer | layout | notes |
|---|---|---|
| `Conv1d` | planar `(C, L)` | zero padding, stride, optional bias; gathers each receptive field then one `dot` per output channel |
| `BatchNorm1d` | planar, in place | folded at load into `scale = w / sqrt(var + eps)`, `shift = b - mean * scale` |
| `Linear` | interleaved `(T, D)` | `gemv` per step |
| `LayerNorm` | interleaved, in place | eps 1e-5 |
| `BiLstm` | interleaved | PyTorch gate order i, f, g, o; `b_ih + b_hh` summed at load; output `[forward, backward]` per step |
| `SelfAttention` | interleaved | `nn.MultiheadAttention` with `batch_first=True`: packed `in_proj`, scaled dot product per head, softmax shifted by its max, `out_proj` |
| `gelu`, `relu`, `sigmoid` | any | GELU is the exact `erf` form, PyTorch's default |

Scratch buffers live in the layers and are sized once, so steady-state
inference allocates nothing. A model instance therefore belongs to one thread.

## Detector

`archive_pipeline.products.detector.WaveformDetector`, `branch1d="cnn-lstm"`,
142k parameters, three seeds averaged in probability.

```
(600, 3) standardised Z/N/E -> asinh
conv 3->24 k7 s2 · BN · GELU      600 -> 300
conv 24->48 k5 s2 · BN · GELU     300 -> 150
conv 48->96 k5 s2 · BN · GELU     150 -> 75
BiLSTM 96 -> 2x48
self-attention (4 heads) + residual -> LayerNorm -> mean over 75 steps
Linear 96->96, times the learned scalar w1
LayerNorm -> Linear 96->96 -> GELU -> Linear 96->1 -> sigmoid
```

`w1` is kept although nothing is fused: training settled it near, not at,
1.0, and dropping it shifts every logit.

Agreement with PyTorch on six real windows, worst case:

| stage | max abs error |
|---|---:|
| conv stack | 3.8e-6 |
| BiLSTM | 1.3e-6 |
| attention | 9.5e-7 |
| pooled embedding | 6.0e-7 |
| logit | 6.0e-7 |
| ensemble probability | identical to 5 decimals |

## Picker

sphase `PhasePicker(arm="wave")`: a dense `{noise, P, S}` distribution over
250 chunks of a 60 s window, so each chunk is 0.24 s.

```
(3, 6000) z-scored
conv 3->32 k9 s2 · BN · ReLU      6000 -> 3000
conv 32->64 k9 s2 · BN · ReLU     3000 -> 1500
conv 64->64 k7 s2 · BN · ReLU     1500 -> 750
average pool 750 -> 250 (exactly 3 to 1)
BiLSTM 64 -> 2x32 -> self-attention + residual -> LayerNorm
Linear 64->64 -> ReLU -> Linear 64->3, per chunk
P, S = argmax chunk of each class, at the chunk centre
```

Logits agree with PyTorch to 5.7e-6, and the picks are identical. On the
M4.9 at DEMI it places P and S 6.96 s apart, which the S-P rule turns into about
58 km. The catalogue distance is 55 km.

**The picker also "picks" in pure noise**, with probabilities around 0.7. It is
therefore only ever run after the detector has triggered.
