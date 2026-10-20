#!/usr/bin/env python3
"""
Plot the results of ./legendre_stability (CSV files in output/) into plots/.

  plots/legendre_P_relative_error.png   |relative error| vs l of P_l(x), |x| < 1  (homework point 5)
  plots/legendre_Q_relative_error.png   |relative error| vs l of Q_l(x),  x > 1
  plots/legendre_values.png             the values themselves: forward, Miller, reference

Usage: python3 plot_legendre.py   (after running ./legendre_stability)
"""

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTPUT_DIR = Path("output")
PLOT_DIR = Path("plots")
EPS = 2.220446049250313e-16  # double precision machine epsilon
FLOOR = 1e-18                # errors that are exactly 0 are drawn at this level

FORWARD = dict(color="tab:blue", marker="o", markersize=3, linewidth=1.5, label="forward")
BACKWARD = dict(color="tab:red", marker="s", markersize=3, linewidth=1.5, label="backward (Miller)")
REFERENCE = dict(color="gray", linewidth=5, alpha=0.5, label="reference (quad)")


def read(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    data = {key: [float(row[key]) for row in rows] for key in rows[0]}
    data["l"] = [int(v) for v in data["l"]]
    return data


def abs_or_floor(values):
    """|v|, with 0 replaced by FLOOR (log scale) and NaN kept as a gap."""
    return [max(abs(v), FLOOR) if not math.isnan(v) else math.nan for v in values]


def grid(n):
    ncols = 2 if n > 1 else 1
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 4.2 * nrows),
                             squeeze=False, constrained_layout=True)
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    return fig, list(axes.flat)


def error_figure(files, function, filename):
    if not files:
        return
    fig, axes = grid(len(files))
    fig.suptitle(f"Relative error of ${function}_\\ell(x)$ in double precision "
                 "(reference: quad precision)", fontsize=14)
    for ax, path in zip(axes, files):
        d = read(path)
        ax.semilogy(d["l"], abs_or_floor(d["rel_err_forward"]), **FORWARD)
        ax.semilogy(d["l"], abs_or_floor(d["rel_err_backward"]), **BACKWARD)
        if function == "P":
            # Error that rounding x itself would cause: kappa * eps (conditioning, not the algorithm)
            ax.semilogy(d["l"], abs_or_floor([c * EPS for c in d["condition_number"]]),
                        color="gray", linestyle="--", linewidth=1.2,
                        label=r"$\kappa_\ell(x)\,\varepsilon$ (conditioning)")
        ax.axhline(EPS, color="black", linestyle=":", linewidth=1, label=r"$\varepsilon$ = 2.2e-16")
        ax.axhline(1.0, color="black", linewidth=0.8, alpha=0.5)
        ax.set_title(f"x = {d['x'][0]:g}")
        ax.set_xlabel(r"degree $\ell$")
        ax.set_ylabel("|relative error|")
        ax.set_ylim(FLOOR / 3, None)
        ax.grid(True, which="major", alpha=0.3)
        ax.legend(fontsize=8, loc="best")
    fig.savefig(PLOT_DIR / filename, dpi=130)
    plt.close(fig)
    print(f"wrote {PLOT_DIR / filename}")


def values_figure(files, filename):
    fig, axes = grid(len(files))
    fig.suptitle("Values: forward and Miller recurrences vs the reference", fontsize=14)
    for ax, path in zip(axes, files):
        d = read(path)
        function = "P" if "_P_" in path.name else "Q"
        if function == "P":
            ax.plot(d["l"], d["reference_quad"], **REFERENCE)
            ax.plot(d["l"], d["forward"], **FORWARD)
            ax.plot(d["l"], d["backward_miller"], **BACKWARD)
            ax.set_title(f"$P_\\ell(x)$, x = {d['x'][0]:g}")
        else:
            # Q_l decreases by many orders of magnitude: plot |Q_l| on a log scale
            # (the wrong forward values also change sign)
            for key, style in (("reference_quad", REFERENCE), ("forward", FORWARD),
                               ("backward_miller", BACKWARD)):
                ax.semilogy(d["l"], [abs(v) if v != 0 else math.nan for v in d[key]], **style)
            ax.set_title(f"$|Q_\\ell(x)|$, x = {d['x'][0]:g}")
        ax.set_xlabel(r"degree $\ell$")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="best")
    fig.savefig(PLOT_DIR / filename, dpi=130)
    plt.close(fig)
    print(f"wrote {PLOT_DIR / filename}")


def by_x(paths):
    return sorted(paths, key=lambda p: read(p)["x"][0])


def main():
    p_files = by_x(OUTPUT_DIR.glob("legendre_P_x*.csv"))
    q_files = by_x(OUTPUT_DIR.glob("legendre_Q_x*.csv"))
    if not p_files and not q_files:
        raise SystemExit("No CSV files in output/: run ./legendre_stability first.")
    PLOT_DIR.mkdir(exist_ok=True)
    error_figure(p_files, "P", "legendre_P_relative_error.png")
    error_figure(q_files, "Q", "legendre_Q_relative_error.png")
    values_figure(p_files + q_files, "legendre_values.png")


if __name__ == "__main__":
    main()
