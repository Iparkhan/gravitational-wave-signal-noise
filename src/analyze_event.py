"""Download and analyze the public 32-second GW150914 H1/L1 strain files."""

from __future__ import annotations

from pathlib import Path
import urllib.request

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import butter, correlate, sosfiltfilt, spectrogram, welch


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
FIGURES = ROOT / "reports" / "figures"
TRIGGER_GPS = 1126259462.422
FILES = {
    "H1": "https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/H-H1_GWOSC_4KHZ_R1-1126259447-32.hdf5",
    "L1": "https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/L-L1_GWOSC_4KHZ_R1-1126259447-32.hdf5",
}
TEMPLATE_URL = "https://gwosc.org/s/events/GW150914/GW150914_4_NR_waveform.txt"


def fetch_files() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for detector, url in FILES.items():
        target = RAW / url.rsplit("/", 1)[-1]
        if not target.exists():
            print(f"Downloading {detector}: {target.name}")
            urllib.request.urlretrieve(url, target)
    template_path = RAW / "GW150914_4_NR_waveform.txt"
    if not template_path.exists():
        print("Downloading numerical-relativity waveform template")
        urllib.request.urlretrieve(TEMPLATE_URL, template_path)


def _decode_strings(values: np.ndarray) -> list[str]:
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def _all_flags_pass(group: h5py.Group, mask_name: str, names_name: str) -> dict[str, bool]:
    names = _decode_strings(np.asarray(group[names_name]))
    masks = np.asarray(group[mask_name], dtype=np.uint32)
    return {name: bool(np.all((masks & (1 << bit)) != 0)) for bit, name in enumerate(names)}


def load_strain(path: Path) -> tuple[np.ndarray, np.ndarray, float, dict[str, bool]]:
    """Return GPS times, strain, sample rate, and segment-wide quality checks."""
    with h5py.File(path, "r") as h5:
        dataset = h5["strain/Strain"]
        strain = np.asarray(dataset, dtype=float)
        gps_start = float(dataset.attrs["Xstart"])
        sample_rate = 1.0 / float(dataset.attrs["Xspacing"])
        dq = _all_flags_pass(h5["quality/simple"], "DQmask", "DQShortnames")
        injections = _all_flags_pass(h5["quality/injections"], "Injmask", "InjShortnames")
        quality = {
            "data_quality_required_flags_pass": all(dq.get(key, False) for key in ("DATA", "CBC_CAT1", "CBC_CAT2")),
            "no_hardware_injection_flags_pass": all(
                value for key, value in injections.items() if key.startswith("NO_")
            ),
            "injection_flags_not_clear": ", ".join(
                key for key, value in injections.items() if key.startswith("NO_") and not value
            ) or "none",
        }
    times = gps_start + np.arange(strain.size) / sample_rate
    return times, strain, sample_rate, quality


def estimate_psd(strain: np.ndarray, sample_rate: float) -> tuple[np.ndarray, np.ndarray]:
    """Estimate the one-sided PSD using one-second Welch windows."""
    return welch(strain, fs=sample_rate, nperseg=int(sample_rate),
                 noverlap=int(sample_rate // 2))


def whiten(strain: np.ndarray, sample_rate: float, psd_freqs: np.ndarray,
            psd: np.ndarray) -> np.ndarray:
    """Whiten strain using an interpolated one-sided PSD, following GWOSC."""
    data_freqs = np.fft.rfftfreq(strain.size, d=1.0 / sample_rate)
    interpolated_psd = np.interp(data_freqs, psd_freqs, psd)
    scale = np.sqrt(np.maximum(interpolated_psd, np.finfo(float).tiny) * sample_rate / 2.0)
    whitened_freq = np.fft.rfft(strain) / scale
    return np.fft.irfft(whitened_freq, n=strain.size)


def main() -> None:
    fetch_files()
    FIGURES.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, float | str]] = []
    traces: dict[str, tuple[np.ndarray, np.ndarray, float]] = {}
    background_distributions: dict[str, tuple[np.ndarray, np.ndarray, float, float]] = {}
    template_time, template_strain = np.loadtxt(RAW / "GW150914_4_NR_waveform.txt").T
    template_zero = int(np.argmin(np.abs(template_time)))
    template_matches: dict[str, tuple[np.ndarray, float, float, int, np.ndarray, np.ndarray]] = {}

    for detector, url in FILES.items():
        path = RAW / url.rsplit("/", 1)[-1]
        times, strain, sample_rate, quality = load_strain(path)
        offset = times - TRIGGER_GPS
        sos = butter(4, [20, 300], btype="bandpass", fs=sample_rate, output="sos")
        psd_freqs, psd = estimate_psd(strain, sample_rate)
        whitened_filtered = sosfiltfilt(sos, whiten(strain, sample_rate, psd_freqs, psd))
        raw_sos = butter(4, [35, 350], btype="bandpass", fs=sample_rate, output="sos")
        filtered = sosfiltfilt(raw_sos, strain)
        window_samples = int(round(0.20 * sample_rate))
        event_center = int(np.argmin(np.abs(offset)))
        event_start = event_center - window_samples // 2
        event_window = filtered[event_start:event_start + window_samples]
        event_rms = float(np.sqrt(np.mean(event_window ** 2)))
        event_peak = float(np.max(np.abs(event_window)))

        step_samples = int(round(0.10 * sample_rate))
        starts = np.arange(0, filtered.size - window_samples + 1, step_samples)
        centers = starts + window_samples // 2
        background_starts = starts[np.abs(offset[centers]) >= 1.0]
        background_windows = np.stack(
            [filtered[start:start + window_samples] for start in background_starts]
        )
        background_rms = np.sqrt(np.mean(background_windows ** 2, axis=1))
        background_peaks = np.max(np.abs(background_windows), axis=1)
        rms_percentile = float(100 * np.mean(background_rms <= event_rms))
        peak_percentile = float(100 * np.mean(background_peaks <= event_peak))

        white_event = whitened_filtered[event_start:event_start + window_samples]
        white_background = np.stack(
            [whitened_filtered[start:start + window_samples] for start in background_starts]
        )
        white_rms = float(np.sqrt(np.mean(white_event ** 2)))
        white_background_rms = np.sqrt(np.mean(white_background ** 2, axis=1))
        white_peak = float(np.max(np.abs(white_event)))
        white_background_peaks = np.max(np.abs(white_background), axis=1)
        background_distributions[detector] = (
            white_background_rms, white_background_peaks, white_rms, white_peak
        )
        rows.append({
            "detector": detector,
            "sample_rate_hz": sample_rate,
            "samples": int(strain.size),
            **quality,
            "event_window_rms_strain": event_rms,
            "background_median_rms_strain": float(np.median(background_rms)),
            "event_rms_percentile_vs_background": rms_percentile,
            "event_peak_abs_strain": event_peak,
            "background_95th_peak_abs_strain": float(np.percentile(background_peaks, 95)),
            "event_peak_percentile_vs_background": peak_percentile,
            "whitened_event_rms_sigma": white_rms,
            "whitened_background_median_rms_sigma": float(np.median(white_background_rms)),
            "whitened_rms_percentile_vs_background": float(100 * np.mean(white_background_rms <= white_rms)),
            "whitened_event_peak_sigma": white_peak,
            "whitened_peak_percentile_vs_background": float(100 * np.mean(white_background_peaks <= white_peak)),
            "background_windows": int(background_windows.shape[0]),
        })
        traces[detector] = (offset, whitened_filtered, sample_rate)

        # Compare the whitened detector data to the public numerical-relativity waveform.
        template_white = sosfiltfilt(sos, whiten(template_strain, sample_rate, psd_freqs, psd))
        match_data = whitened_filtered.copy()
        if detector == "L1":
            match_data = -np.roll(match_data, int(round(0.007 * sample_rate)))
        template_len = template_white.size
        if template_len > match_data.size:
            raise ValueError("Waveform template is longer than the detector data segment")
        raw_corr = correlate(match_data, template_white, mode="valid", method="fft")
        cumulative_energy = np.concatenate(([0.0], np.cumsum(match_data ** 2)))
        local_energy = cumulative_energy[template_len:] - cumulative_energy[:-template_len]
        template_norm = float(np.linalg.norm(template_white))
        normalized_corr = raw_corr / np.maximum(template_norm * np.sqrt(local_energy), 1e-30)
        all_starts = np.arange(normalized_corr.size)
        merge_offsets = offset[all_starts + template_zero]
        match_step = int(round(0.10 * sample_rate))
        event_lags = np.flatnonzero(np.abs(merge_offsets) <= 0.05)
        background_lags = np.arange(0, normalized_corr.size, match_step)
        background_lags = background_lags[np.abs(merge_offsets[background_lags]) > 1.0]
        event_lag = int(event_lags[np.argmax(np.abs(normalized_corr[event_lags]))])
        background_scores = np.abs(normalized_corr[background_lags])
        event_score = float(abs(normalized_corr[event_lag]))
        match_percentile = float(100 * np.mean(background_scores <= event_score))
        rows[-1].update({
            "template_match_abs_correlation": event_score,
            "template_match_event_offset_s": float(merge_offsets[event_lag]),
            "template_match_percentile_vs_background": match_percentile,
            "template_background_windows": int(background_scores.size),
        })
        template_segment = template_white
        data_segment = match_data[event_lag:event_lag + template_len]
        best_scale = float(np.dot(data_segment, template_segment) /
                           max(np.dot(template_segment, template_segment), 1e-30))
        template_offset = template_time - template_time[template_zero] + merge_offsets[event_lag]
        data_offset = times[event_lag:event_lag + template_len] - TRIGGER_GPS
        template_matches[detector] = (
            merge_offsets, normalized_corr, event_score, event_lag,
            data_segment, best_scale * template_segment
        )

        fig, ax = plt.subplots(figsize=(10, 4.5))
        focus = (data_offset >= -0.15) & (data_offset <= 0.05)
        ax.plot(data_offset[focus], data_segment[focus], lw=0.9, label=f"{detector} whitened data")
        template_focus = (template_offset >= -0.15) & (template_offset <= 0.05)
        ax.plot(template_offset[template_focus], (best_scale * template_segment)[template_focus],
                lw=1.2, label="Scaled numerical-relativity template")
        ax.axvline(0, color="black", ls="--", lw=1, label="Published event time")
        ax.set(title=f"Waveform template comparison — {detector}",
               xlabel="Time from published GPS trigger (s)", ylabel="Whitened, band-passed strain (σ)")
        ax.legend()
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(FIGURES / f"{detector.lower()}_template_overlay.png", dpi=160)
        plt.close(fig)

        near = np.abs(offset) <= 0.25
        fig, ax = plt.subplots(figsize=(10, 4.5))
        ax.plot(offset[near], whitened_filtered[near], lw=0.8)
        ax.axvline(0, color="crimson", ls="--", lw=1, label="Published event time")
        ax.set(title=f"GW150914 around event time — {detector}",
               xlabel="Time from published GPS trigger (s)", ylabel="Whitened, band-passed strain (σ)")
        ax.legend()
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(FIGURES / f"{detector.lower()}_event_window.png", dpi=160)
        plt.close(fig)

        f, spec_t, power = spectrogram(filtered, fs=sample_rate, nperseg=256,
                                       noverlap=224, scaling="density", mode="psd")
        # Normalize each frequency bin by its median level in the full 32 s file
        # so persistent detector lines do not dominate the event-centered view.
        relative_power = power / np.maximum(np.median(power, axis=1, keepdims=True), 1e-60)
        time_mask = np.abs(spec_t + times[0] - TRIGGER_GPS) <= 0.5
        freq_mask = (f >= 20) & (f <= 500)
        fig, ax = plt.subplots(figsize=(10, 4.5))
        mesh = ax.pcolormesh(spec_t[time_mask] + times[0] - TRIGGER_GPS,
                             f[freq_mask], 10 * np.log10(relative_power[freq_mask][:, time_mask] + 1e-60),
                             shading="auto", cmap="magma")
        ax.axvline(0, color="cyan", ls="--", lw=1)
        ax.set(title=f"GW150914 time-frequency view — {detector}",
               xlabel="Time from published GPS trigger (s)", ylabel="Frequency (Hz)")
        fig.colorbar(mesh, ax=ax, label="Power relative to each bin's 32 s median (dB)")
        fig.tight_layout()
        fig.savefig(FIGURES / f"{detector.lower()}_spectrogram.png", dpi=160)
        plt.close(fig)

    # Compare the event's whitened features with the empirical off-event distributions.
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    detectors = ["H1", "L1"]
    colors = ["#2878b5", "#e87522"]
    rms_samples = [background_distributions[d][0] for d in detectors]
    peak_samples = [background_distributions[d][1] for d in detectors]
    axes[0].boxplot(rms_samples, tick_labels=detectors, showfliers=False)
    axes[1].boxplot(peak_samples, tick_labels=detectors, showfliers=False)
    for position, (detector, color) in enumerate(zip(detectors, colors), start=1):
        _, _, event_rms, event_peak = background_distributions[detector]
        axes[0].scatter(position, event_rms, color=color, marker="D", s=55, zorder=3,
                        label=f"{detector} event ({event_rms:.2f}σ)")
        axes[1].scatter(position, event_peak, color=color, marker="D", s=55, zorder=3,
                        label=f"{detector} event ({event_peak:.2f}σ)")
    axes[0].set(title="Equal-window whitened RMS", xlabel="Detector", ylabel="Window RMS (σ)")
    axes[1].set(title="Equal-window whitened peak", xlabel="Detector", ylabel="Maximum absolute value (σ)")
    for ax in axes:
        ax.set_yscale("log")
        ax.grid(alpha=0.2, which="both")
        ax.legend(fontsize=8)
    fig.suptitle("Event window compared with off-event windows (298 per detector)")
    fig.tight_layout()
    fig.savefig(FIGURES / "background_comparison.png", dpi=160)
    plt.close(fig)

    # The GWOSC educational tutorial aligns this event by shifting L1 +7 ms
    # and inverting its sign. This is a visual comparison, not a new delay fit.
    h1_offset, h1_white, sample_rate = traces["H1"]
    l1_offset, l1_white, _ = traces["L1"]
    shift_samples = int(round(0.007 * sample_rate))
    l1_aligned = -np.roll(l1_white, shift_samples)
    near = (h1_offset >= -0.15) & (h1_offset <= 0.05)
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(h1_offset[near], h1_white[near], lw=0.9, label="H1")
    ax.plot(l1_offset[near], l1_aligned[near], lw=0.9, label="L1 (shifted +7 ms, sign inverted)")
    ax.axvline(0, color="black", ls="--", lw=1, label="Published event time")
    ax.set(title="Whitened GW150914 strain in two detectors",
           xlabel="Time from published GPS trigger (s)", ylabel="Whitened, band-passed strain (σ)")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIGURES / "detector_coincidence.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 4.5))
    for detector, (_, normalized_corr, event_score, event_lag, _, _) in template_matches.items():
        offset, _, _ = traces[detector]
        match_offsets = offset[np.arange(normalized_corr.size) + template_zero]
        focus = np.abs(match_offsets) <= 0.15
        ax.plot(match_offsets[focus], np.abs(normalized_corr[focus]), lw=0.9,
                label=f"{detector} |r| (event max {event_score:.2f})")
        ax.scatter(match_offsets[event_lag], event_score, s=35)
    ax.axvline(0, color="black", ls="--", lw=1, label="Published event time")
    ax.set(title="Normalized waveform-template correlation near GW150914",
           xlabel="Template merger time relative to published trigger (s)",
           ylabel="Absolute normalized correlation")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIGURES / "template_match_near_event.png", dpi=160)
    plt.close(fig)

    summary = pd.DataFrame(rows)
    summary.to_csv(ROOT / "reports" / "detector_summary.csv", index=False)
    print(summary.to_string(index=False))
    print(f"\nFigures saved in {FIGURES}")
    print("Interpretation caution: this compares one known event window with nearby data; it is not a detection-significance test.")


if __name__ == "__main__":
    main()
