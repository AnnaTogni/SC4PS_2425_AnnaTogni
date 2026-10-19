#!/usr/bin/env python3

"""
Plot the coupled-oscillator data with matplotlib.

This script complements the C simulation by turning the generated CSV files into
figures suitable for inspection or discussion. It also overlays the theoretical
normal-mode frequencies on the FFT figure, which is a useful research habit:
always compare numerical output with theoretical expectations as directly as
possible.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / "output"
PLOT_DIR = ROOT / "plots"

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


def read_two_series(path: Path) -> tuple[list[float], list[float], list[float]]:
    """Read a CSV file containing one x column and two y columns."""
    x_values: list[float] = []
    y1_values: list[float] = []
    y2_values: list[float] = []

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            x_values.append(float(row[0]))
            y1_values.append(float(row[1]))
            y2_values.append(float(row[2]))

    return x_values, y1_values, y2_values


def read_csv_columns(path: Path) -> tuple[list[float], list[float]]:
    """Read a two-column CSV file after skipping its header row."""
    x_values: list[float] = []
    y_values: list[float] = []

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            x_values.append(float(row[0]))
            y_values.append(float(row[1]))

    return x_values, y_values


def plot_time_series() -> None:
    """Plot the two displacements in time and a zoomed early-time view."""
    time_s, x1_m, x2_m = read_two_series(OUTPUT_DIR / "coupled_oscillators_time.csv")

    fig, axes = plt.subplots(2, 1, figsize=(14, 7), constrained_layout=True)
    fig.suptitle("Coupled oscillators: displacement versus time")

    axes[0].plot(time_s, x1_m, linewidth=1.8, label="x1(t)")
    axes[0].plot(time_s, x2_m, linewidth=1.8, label="x2(t)")
    axes[0].set_title("Mass 1 and mass 2")
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Displacement (m)")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    early_time = [index for index, value in enumerate(time_s) if value <= 20.0]
    last_index = early_time[-1] + 1 if early_time else len(time_s)
    axes[1].plot(time_s[:last_index], x1_m[:last_index], linewidth=1.8, label="x1(t)")
    axes[1].set_title("Mass 1 over an early-time window")
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Displacement (m)")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    fig.savefig(PLOT_DIR / "coupled_oscillators_time.png", dpi=150)
    plt.close(fig)


def plot_spectrum() -> None:
    """Plot the FFT spectrum together with theoretical normal-mode markers."""
    frequency_hz, amplitude = read_csv_columns(
        OUTPUT_DIR / "coupled_oscillators_spectrum.csv"
    )
    mass = 1.0
    wall_k = 25.0
    coupling_k = 7.0
    in_phase_hz = math.sqrt(wall_k / mass) / (2.0 * math.pi)
    out_of_phase_hz = math.sqrt((wall_k + 2.0 * coupling_k) / mass) / (2.0 * math.pi)

    fig, axis = plt.subplots(figsize=(14, 5), constrained_layout=True)
    axis.plot(frequency_hz, amplitude, linewidth=1.8, label="FFT amplitude")
    axis.axvline(
        in_phase_hz,
        color="tab:red",
        linestyle="--",
        linewidth=1.8,
        label=f"Theory in-phase: {in_phase_hz:.4f} Hz",
    )
    axis.axvline(
        out_of_phase_hz,
        color="tab:green",
        linestyle="--",
        linewidth=1.8,
        label=f"Theory out-of-phase: {out_of_phase_hz:.4f} Hz",
    )
    axis.set_xlim(0.0, 2.0)
    axis.set_title("Coupled oscillators: FFT amplitude spectrum of x1(t)")
    axis.set_xlabel("Frequency (Hz)")
    axis.set_ylabel("Amplitude")
    axis.grid(True, alpha=0.3)
    axis.legend()

    fig.savefig(PLOT_DIR / "coupled_oscillators_spectrum.png", dpi=150)
    plt.close(fig)


def main() -> None:
    """Create the output directory and generate both coupled-oscillator figures."""
    PLOT_DIR.mkdir(exist_ok=True)
    plot_time_series()
    plot_spectrum()


if __name__ == "__main__":
    main()
