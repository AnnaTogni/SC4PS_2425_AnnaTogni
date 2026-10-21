#!/usr/bin/env python3
"""
Figures for HomeWork7 from the CSV files in output/ (run the five programs first).

  plots/coin_tosses.png      running fraction of heads with the +-1/(2 sqrt(n)) band
  plots/monte_carlo_pi.png   |error| of the pi estimate vs N, with the 1/sqrt(N) law
  plots/square_uniform.png   histogram of Y = U^2 vs f(y) = 1/(2 sqrt(y))
  plots/exponential.png      histogram vs lambda exp(-lambda y), linear and log scale
  plots/empirical_cdf.png    empirical vs exact CDF, their difference, and D_n vs n
  plots/randu_planes.png     bonus: RANDU's triples on 15 planes vs the notebook LCG

Usage: python3 plot_random.py
"""

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("output")
PLOTS = Path("plots")
LAMBDA = 1.5


def read(name):
    with (OUT / name).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {key: [float(r[key]) for r in rows] for key in rows[0]}


def save(fig, name):
    fig.savefig(PLOTS / name, dpi=130)
    plt.close(fig)
    print(f"wrote {PLOTS / name}")


def coin_tosses():
    d = read("coin_tosses.csv")
    n, f, s = d["n"], d["fraction_heads"], d["sigma"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)
    fig.suptitle("Exercise 1: running fraction of heads (law of large numbers)")
    for ax, log in zip(axes, (True, False)):
        ax.plot(n, f, color="tab:blue", linewidth=1, label="running fraction $f_n$")
        # +-2 sigma band, clipped to [0, 1] (a fraction cannot leave that interval)
        ax.fill_between(n, [max(0.0, 0.5 - 2 * v) for v in s], [min(1.0, 0.5 + 2 * v) for v in s],
                        color="tab:orange", alpha=0.25, label=r"$0.5 \pm 2\sigma_n$, $\sigma_n = 1/(2\sqrt{n})$")
        ax.axhline(0.5, color="black", linewidth=1)
        ax.set_xlabel("number of tosses n")
        ax.set_ylabel("fraction of heads")
        ax.grid(True, alpha=0.3)
        if log:
            ax.set_xscale("log")
            ax.set_title("log scale: whole run")
        else:
            ax.set_ylim(0.48, 0.52)
            ax.set_title("linear scale, zoom around 0.5")
        ax.legend(loc="upper right", fontsize=9)
    save(fig, "coin_tosses.png")


def monte_carlo_pi():
    d = read("monte_carlo_pi.csv")
    N = d["N"]
    fig, ax = plt.subplots(figsize=(7.5, 5), constrained_layout=True)
    ax.loglog(N, [abs(e) for e in d["error_first_run"]], "o", color="tab:blue",
              label=r"$|\pi_N - \pi|$, single run")
    ax.loglog(N, d["rms_error"], "s-", color="tab:red", label="RMS error over repeated runs")
    ax.loglog(N, d["expected_sigma"], "--", color="black",
              label=r"theory $\sqrt{\pi(4-\pi)/N} \approx 1.64/\sqrt{N}$")
    ax.set_xlabel("number of points N")
    ax.set_ylabel(r"error of the estimate of $\pi$")
    ax.set_title("Exercise 2: Monte Carlo estimate of $\\pi$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    save(fig, "monte_carlo_pi.png")


def histogram_figure(name, title, pdf_label, pdf, log_panel):
    d = read(name)
    left, right, dens = d["bin_left"], d["bin_right"], d["density_histogram"]
    width = right[0] - left[0]
    ymax = right[-1]
    fine = [ymax * (i + 0.5) / 2000 for i in range(2000)]
    ncols = 2 if log_panel else 1
    fig, axes = plt.subplots(1, ncols, figsize=(6.5 * ncols, 4.8), constrained_layout=True,
                             squeeze=False)
    fig.suptitle(title)
    for k, ax in enumerate(axes[0]):
        ax.bar(left, dens, width=width, align="edge", color="tab:blue", alpha=0.45,
               edgecolor="tab:blue", label="histogram (normalized)")
        ax.plot(fine, [pdf(y) for y in fine], color="tab:red", linewidth=2, label=pdf_label)
        ax.step(left + [right[-1]], d["density_bin_average"] + [d["density_bin_average"][-1]],
                where="post", color="black", linewidth=1, linestyle=":",
                label="pdf averaged over each bin")
        ax.set_xlabel("y")
        ax.set_ylabel("probability density")
        ax.grid(True, alpha=0.3)
        if k == 1:
            ax.set_yscale("log")
            ax.set_title("log scale")
        ax.legend(fontsize=9)
    return fig, axes[0]


def square_uniform():
    fig, axes = histogram_figure("square_uniform.csv", "Exercise 3: $Y = U^2$, $U \\sim$ Uniform(0,1)",
                                 r"$f_Y(y) = 1/(2\sqrt{y})$", lambda y: 1 / (2 * math.sqrt(y)),
                                 log_panel=True)
    axes[0].set_ylim(0, 8)
    axes[1].set_xscale("log")
    axes[1].set_title("log-log scale: $f_Y \\propto y^{-1/2}$ is a straight line")
    save(fig, "square_uniform.png")


def exponential():
    fig, axes = histogram_figure("exponential.csv",
                                 f"Exercise 4: exponential by inverse transform, $\\lambda = {LAMBDA}$",
                                 r"$f(y) = \lambda e^{-\lambda y}$",
                                 lambda y: LAMBDA * math.exp(-LAMBDA * y), log_panel=True)
    axes[1].set_title("log scale: the exponential is a straight line")
    save(fig, "exponential.png")


def empirical_cdf():
    d = read("empirical_cdf.csv")
    ks = read("ks_vs_n.csv")
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8), constrained_layout=True)
    fig.suptitle(f"Exercise 5: empirical CDF of the exponential sample ($\\lambda = {LAMBDA}$)")

    ax = axes[0]
    ax.step(d["y"], d["ecdf"], where="post", color="tab:blue", linewidth=2.5, label="empirical $F_n(y)$")
    ax.plot(d["y"], d["exact_cdf"], color="tab:red", linewidth=1.2, linestyle="--",
            label=r"exact $F(y) = 1 - e^{-\lambda y}$")
    ax.set_xlabel("y")
    ax.set_ylabel("CDF")
    ax.set_xlim(0, 5)
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_title("the two curves overlap")

    ax = axes[1]
    ax.plot(d["y"], [e - x for e, x in zip(d["ecdf"], d["exact_cdf"])], color="tab:blue", linewidth=1)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("y")
    ax.set_ylabel("$F_n(y) - F(y)$")
    ax.set_xlim(0, 5)
    ax.grid(True, alpha=0.3)
    ax.set_title("difference: its largest absolute value is $D_n$")

    ax = axes[2]
    ax.loglog(ks["n"], ks["D_n"], "o-", color="tab:blue", label=r"$D_n = \max_y |F_n - F|$")
    ax.loglog(ks["n"], [1.36 / math.sqrt(m) for m in ks["n"]], "--", color="tab:red",
              label=r"95% limit $1.36/\sqrt{n}$")
    ax.set_xlabel("sample size n")
    ax.set_ylabel("Kolmogorov-Smirnov distance")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    ax.set_title("$D_n$ shrinks like $1/\\sqrt{n}$")
    save(fig, "empirical_cdf.png")


def randu_planes():
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers the 3D projection)

    r, l, comb = read("randu_triples.csv"), read("lcg_triples.csv"), read("plane_combination.csv")
    fig = plt.figure(figsize=(17, 5.2), constrained_layout=True)
    fig.suptitle("Bonus: consecutive triples $(u_n, u_{n+1}, u_{n+2})$ -- RANDU vs notebook LCG")
    for k, (d, name, color) in enumerate(((r, "RANDU", "tab:red"), (l, "notebook LCG", "tab:blue"))):
        ax = fig.add_subplot(1, 3, k + 1, projection="3d")
        ax.scatter(d["u0"], d["u1"], d["u2"], s=2, alpha=0.5, color=color)
        # look along (2, 3, 0), a direction lying in RANDU's planes 9x - 6y + z = k:
        # the planes are then seen edge-on, as lines
        ax.view_init(elev=0, azim=math.degrees(math.atan2(3, 2)))
        ax.set_title(f"{name} (3000 triples, seen along the planes)")
        ax.set_xlabel("$u_n$")
        ax.set_ylabel("$u_{n+1}$")
        ax.set_zlabel("$u_{n+2}$")
    ax = fig.add_subplot(1, 3, 3)
    ax.hist(comb["s_lcg"], bins=400, range=(-6, 10), color="tab:blue", alpha=0.6, label="notebook LCG")
    ax.hist(comb["s_randu"], bins=400, range=(-6, 10), color="tab:red", alpha=0.8, label="RANDU")
    ax.set_yscale("log")
    ax.set_xlabel(r"$s = 9u_n - 6u_{n+1} + u_{n+2}$")
    ax.set_ylabel("count (100000 triples)")
    ax.set_title("RANDU: s is always an integer (15 values)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    save(fig, "randu_planes.png")


def main():
    PLOTS.mkdir(exist_ok=True)
    coin_tosses()
    monte_carlo_pi()
    square_uniform()
    exponential()
    empirical_cdf()
    randu_planes()


if __name__ == "__main__":
    main()
