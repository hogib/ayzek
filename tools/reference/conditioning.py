"""Reference conditioning: what the models' training data went through.

The C++ tests compare against these, so they are copies, not rewrites:

- `taper_vector`, `_bandpass`, `clean_block`: the detector and picker windows,
  from archive_pipeline (`archive/clean.py`, at db13814)
- `standardize`, `clean_and_filter_1d`, `SpectrogramEncoder`: the magnitude
  windows, from data_downloader (`seismic_cli/core.py` and `spectrogram.py`, at
  66f3cac), the encoder that built the regressor's training set
"""
import math
from typing import Dict, Optional, Tuple

import numpy as np
from scipy import signal


def taper_vector(n):
    """The 5% Hann cosine taper applied to each end of a window of length `n`."""
    t = np.ones(n)
    k = int(n * 0.05)
    if k > 0:
        w = signal.windows.hann(k * 2)
        t[:k] = w[:k]
        t[-k:] = w[k:]
    return t


def _bandpass(fs, freqmin, freqmax):
    """Butterworth coefficients, or `(None, None)` when the band is degenerate.

    A band whose high corner has been clamped to below its low corner cannot be
    filtered; the samples pass through untouched rather than the call failing,
    which is what happens at a station recorded below 2*freqmax.
    """
    nyquist = fs / 2.0
    high = freqmax if nyquist > freqmax else nyquist - 1.0
    if high <= freqmin:
        return None, None
    return signal.butter(4, [freqmin, high], btype="bandpass", fs=fs)


def clean_block(x, fs, freqmin, freqmax, taper):
    """Conditions an `(n, win)` block of windows at once.

    Identical in operation and order to `clean_window`; scipy is given an axis
    instead of being called once per row, which is the whole reason a block
    form exists.

    Args:
        x: `(n, win)` array of raw samples. Not modified in place.
        fs: Sampling rate in Hz.
        freqmin: Bandpass low corner in Hz.
        freqmax: Bandpass high corner in Hz, clamped below Nyquist.
        taper: Taper of length `win`, from `taper_vector`.

    Returns:
        `(n, win)` array of conditioned samples.
    """
    x = signal.detrend(x, type="linear", axis=-1)
    x = signal.detrend(x, type="constant", axis=-1)
    x = x * taper
    b, a = _bandpass(fs, freqmin, freqmax)
    return x if b is None else signal.filtfilt(b, a, x, axis=-1)


# --- magnitude windows ------------------------------------------------------

def standardize(x: np.ndarray, mu: Optional[float] = None, sigma: Optional[float] = None,
                 eps: float = 1e-12) -> np.ndarray:
    """
    Standardizes x using EITHER a provided (mu, sigma) -- e.g. a station's
    long-term noise baseline -- OR, if not provided, the window's own
    mean/std (plain per-window self-standardization).
    """
    x = np.asarray(x, dtype=np.float64)
    if mu is None:
        mu = np.mean(x)
    if sigma is None:
        sigma = np.std(x)
    if sigma < eps:
        sigma = eps
    return (x - mu) / sigma


def clean_and_filter_1d(x: np.ndarray, fs: float, freqmin: float, freqmax: float) -> np.ndarray:
    x = signal.detrend(x, type='linear')
    x = signal.detrend(x, type='constant')
    n = len(x)
    taper_len = int(n * 0.05)
    if taper_len > 0:
        window = signal.windows.hann(taper_len * 2)
        x[:taper_len] *= window[:taper_len]
        x[-taper_len:] *= window[-taper_len:]

    nyquist = fs / 2.0
    actual_freqmax = freqmax if nyquist > freqmax else nyquist - 1.0

    if actual_freqmax > freqmin:
        b, a = signal.butter(4, [freqmin, actual_freqmax], btype='bandpass', fs=fs)
        x = signal.filtfilt(b, a, x)

    return x


# CHANNEL SELECTION

# Preference order per component role. Sorting channel letters alphabetically
# and taking the first three could grab e.g. ['1','2','E'] at a station with
# mixed sensor codes -- two horizontals from one instrument plus one from
# another, and no vertical at all. Selecting by explicit role keeps the
# component->color mapping fixed (R=Z, G=N-ish, B=E-ish) for every station.
_COMPONENT_ROLES = (('Z',), ('N', '1'), ('E', '2'))


NORMALIZE_MODES = ("station", "per_window", "none")

# An STFT frame spanning more than this fraction of the analysis window has too
# little of the window left to resolve anything in time. 0.25 keeps at least
# four independent frame positions; 0.50 and above is effectively a single
# smeared frame and is flagged harder.
NFFT_WINDOW_FRACTION_WARN = 0.25
NFFT_WINDOW_FRACTION_SEVERE = 0.50
MIN_USEFUL_FRAMES = 8
# Below this, frequency resolution (fs / n_fft) gets too coarse to separate the
# 1-45 Hz band regardless of how many frames it buys.
MIN_SUGGESTED_NFFT = 32


def stft_geometry(window_seconds: float, nominal_fs: float, n_fft: int,
                  hop_length: Optional[int] = None) -> Dict[str, float]:
    """
    Resolves what an (n_fft, hop) choice actually produces for a given window.

    `frames` matches torchaudio's `Spectrogram` with its default `center=True`
    (`n_samples // hop + 1`), verified against emitted tensors at both 3 s and
    6 s rather than assumed.
    """
    hop_length = hop_length if hop_length is not None else n_fft // 4
    n_samples = int(round(nominal_fs * window_seconds))
    return {
        "n_samples": n_samples,
        "hop_length": hop_length,
        "freq_bins": n_fft // 2 + 1,
        "frames": (n_samples // hop_length) + 1 if hop_length > 0 else 0,
        "nfft_seconds": n_fft / nominal_fs,
        "window_fraction": (n_fft / n_samples) if n_samples else float("inf"),
        "freq_resolution_hz": nominal_fs / n_fft,
    }


def suggest_stft_params(window_seconds: float, nominal_fs: float) -> Tuple[int, int]:
    """
    Largest power-of-two `n_fft` spanning at most a quarter of the window, with
    `hop = n_fft // 4`.

    Tying n_fft to the window rather than fixing it means time resolution stays
    constant as the window changes: this yields 19 frames at both 3 s and 6 s,
    where a fixed n_fft=256 gives 10 frames at 6 s and only 5 at 3 s.
    """
    n_samples = int(round(nominal_fs * window_seconds))
    budget = max(MIN_SUGGESTED_NFFT, int(n_samples * NFFT_WINDOW_FRACTION_WARN))
    n_fft = 1 << max(int(budget).bit_length() - 1, 5)   # floor to a power of two, >= 32
    return n_fft, max(1, n_fft // 4)


def check_stft_resolution(window_seconds: float, nominal_fs: float, n_fft: int,
                          hop_length: Optional[int] = None) -> Optional[str]:
    """
    Returns a warning string when (n_fft, hop) is badly matched to the window,
    else None.

    This exists because a 3 s dataset was generated with the 6 s default of
    n_fft=256: one frame then spans 85% of the window and the tensor collapses
    to 5 time bins, which is a silently degraded input rather than an error.
    Nothing in the pipeline noticed.
    """
    g = stft_geometry(window_seconds, nominal_fs, n_fft, hop_length)
    sug_nfft, sug_hop = suggest_stft_params(window_seconds, nominal_fs)
    if g["window_fraction"] <= NFFT_WINDOW_FRACTION_WARN and g["frames"] >= MIN_USEFUL_FRAMES:
        return None
    if n_fft <= sug_nfft:
        # Already at or below what this window can support. On very short
        # windows `MIN_SUGGESTED_NFFT` binds and the recommendation itself
        # exceeds the fraction budget; warning here would demand a value the
        # function cannot recommend, so say nothing rather than something
        # unactionable.
        return None
    severity = "SEVERE" if g["window_fraction"] > NFFT_WINDOW_FRACTION_SEVERE else "WARN"
    sug = stft_geometry(window_seconds, nominal_fs, sug_nfft, sug_hop)
    return (
        f"[{severity}] n_fft={n_fft} is poorly matched to a {window_seconds:g}s window "
        f"at {nominal_fs:g} Hz.\n"
        f"        One frame spans {g['nfft_seconds']:.2f}s = {g['window_fraction'] * 100:.0f}% "
        f"of the window, leaving only {int(g['frames'])} time bins "
        f"(tensor {g['freq_bins']} x {int(g['frames'])}).\n"
        f"        Suggested: --n-fft {sug_nfft} --hop-length {sug_hop}  -> "
        f"{sug['freq_bins']} x {int(sug['frames'])} bins, "
        f"{sug['freq_resolution_hz']:.2f} Hz resolution, "
        f"{sug['window_fraction'] * 100:.0f}% of the window per frame.\n"
        f"        Keeping the current value is fine ONLY if you are matching an "
        f"existing dataset -- a model must see the same geometry it was trained on."
    )


def _resample_to(x: np.ndarray, fs_from: float, fs_to: float) -> np.ndarray:
    """Polyphase resample so every window yields the same number of samples."""
    if abs(fs_from - fs_to) < 1e-9:
        return x
    g = math.gcd(int(round(fs_to)), int(round(fs_from)))
    up, down = int(round(fs_to)) // g, int(round(fs_from)) // g
    return signal.resample_poly(x, up, down, axis=0)


def _fit_length(x: np.ndarray, n: int) -> np.ndarray:
    """Force exactly n samples (resampling can be off by one)."""
    if len(x) == n:
        return x
    if len(x) > n:
        return x[:n]
    return np.pad(x, ((0, n - len(x)),) + ((0, 0),) * (x.ndim - 1), mode="constant")


class SpectrogramEncoder:
    """
    Per-window encoder producing a (3, freq, time) float32 tensor.

    Picklable by design (plain attributes only) so it survives the
    ProcessPoolExecutor hand-off; torch/torchaudio are imported lazily inside
    the worker so the RAM-only path never needs them installed.

    `requires_spawn` tells the orchestrator to use 'spawn' workers. Under the
    default 'fork', torch's threading state does not survive the fork and the
    workers deadlock silently -- sleeping at 0% CPU with nothing written, which
    looks like a slow run rather than a hang. (The standalone script this
    replaces had the same fork+torch structure.)
    """
    ext = ".pt"
    requires_spawn = True

    def __init__(self, n_fft: int = 256, hop_length: Optional[int] = None,
                 top_db: float = 80.0, nominal_fs: float = 100.0,
                 window_seconds: float = 60.0, normalize: str = "station",
                 noise_profiles: Optional[Dict[Tuple[str, str], np.ndarray]] = None):
        if normalize not in NORMALIZE_MODES:
            raise ValueError(f"normalize must be one of {NORMALIZE_MODES}, got {normalize!r}")
        self.n_fft = n_fft
        self.hop_length = hop_length if hop_length is not None else n_fft // 4
        self.top_db = top_db
        self.nominal_fs = nominal_fs
        self.window_seconds = window_seconds
        self.normalize = normalize
        self.noise_profiles = noise_profiles or {}
        self._tf = None  # lazily built per worker process

        # Fires once, where the encoder is constructed -- workers receive it
        # unpickled, so this does not repeat per process.
        problem = check_stft_resolution(window_seconds, nominal_fs, n_fft, hop_length)
        if problem:
            print(problem)

    # -- lazy torch setup -------------------------------------------------
    def _transforms(self):
        if self._tf is None:
            import torch
            import torchaudio.transforms as T
            # One thread per worker: the pool already provides parallelism, and
            # letting each worker spin up a full thread pool oversubscribes the
            # machine badly.
            torch.set_num_threads(1)
            self._tf = (
                torch,
                T.Spectrogram(n_fft=self.n_fft, hop_length=self.hop_length, power=2.0),
                # top_db bounds the dynamic range; without it a near-silent bin
                # drags the tensor's floor arbitrarily low and destabilizes
                # whatever normalization runs afterwards.
                T.AmplitudeToDB(stype="power", top_db=self.top_db),
            )
        return self._tf

    def target_samples(self) -> int:
        return int(round(self.nominal_fs * self.window_seconds))

    def spec_db(self, cleaned_win: np.ndarray, fs_station: float):
        """(samples, 3) float array -> (3, freq, time) dB tensor at nominal fs."""
        torch, spec_tf, db_tf = self._transforms()
        x = _fit_length(_resample_to(cleaned_win, fs_station, self.nominal_fs),
                        self.target_samples())
        t = torch.from_numpy(np.ascontiguousarray(x.T)).float()
        return db_tf(spec_tf(t))

    def normalize_spec(self, s, sta_key, selection):
        """Applies `self.normalize` to a raw (3, freq, time) dB tensor. Split
        out from `__call__` so `SpectrogramDualEncoder` can reuse the exact
        same normalization instead of a second, drifting copy of it."""
        torch, _, _ = self._transforms()
        if self.normalize == "station":
            applied = False
            if self.noise_profiles:
                prof = [self.noise_profiles.get((sta_key, c)) for c in selection]
                if all(p is not None for p in prof):
                    ref = torch.from_numpy(np.stack(prof)).float().unsqueeze(-1)
                    s = s - ref            # dB above this station's noise floor
                    applied = True
            if not applied:               # no usable profile -> safe fallback
                s = (s - s.mean()) / (s.std() + 1e-6)
        elif self.normalize == "per_window":
            s = (s - s.mean()) / (s.std() + 1e-6)
        return s

    # -- encoder protocol -------------------------------------------------
    def __call__(self, cleaned_win, fs_station, sta_key, selection,
                 station_baselines, out_dir, stem):
        torch, _, _ = self._transforms()
        s = self.normalize_spec(self.spec_db(cleaned_win, fs_station), sta_key, selection)
        filename = stem + self.ext
        torch.save(s.contiguous(), out_dir / filename)
        return filename
