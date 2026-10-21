# Random numbers: law of large numbers, Monte Carlo, transformations

Five exercises on pseudo-random numbers, written in C with the GNU Scientific Library.

**The generator.** All programs use GSL's **Mersenne Twister** (`gsl_rng_mt19937`, period
2¹⁹⁹³⁷ − 1) with a fixed seed (default 12345). Every result below is therefore exactly
reproducible, and a different seed gives an independent sample. `gsl_rng_uniform` returns
U ∈ [0, 1).

| File | Content |
|---|---|
| `coin_tosses.c` | 1. Coin tosses and the law of large numbers |
| `monte_carlo_pi.c` | 2. Monte Carlo estimate of π with two RNG streams |
| `square_uniform.c` | 3. Change of variables Y = U² |
| `exponential.c`, `exponential.h` | 4. Exponential distribution by inverse transform |
| `empirical_cdf.c` | 5. Empirical CDF of the sample of exercise 4 |
| `common.c`, `common.h` | Shared helpers: generator, arguments, output files, χ² test |
| `plot_random.py` | Figures in `plots/` |
| `Makefile` | `make` (builds and prints the instructions), `make run`, `make plot`, `make help`, `make clean` |

## Build and run

```bash
make            # builds the five programs and prints how to use them
make run        # runs the five exercises (CSV files in output/)
make plot       # runs them and draws plots/*.png
```

Every program takes optional arguments (sample size etc.); the last one is always the seed.
For example, `./coin_tosses 1e6 7` does 10⁶ tosses with seed 7.

---

## 1. Coin tosses and the law of large numbers

100 000 tosses with `gsl_rng_uniform_int(r, 2)` (0 or 1 with equal probability). Running
fraction of heads fₙ:

| n | Heads | fₙ | \|fₙ − 0.5\| | In units of σₙ |
|---|---|---|---|---|
| 10 | 5 | 0.500000 | 0 | 0.00 |
| 100 | 62 | 0.620000 | 0.12 | 2.40 |
| 1 000 | 500 | 0.500000 | 0 | 0.00 |
| 10 000 | 4 948 | 0.494800 | 0.0052 | 1.04 |
| 100 000 | 49 952 | 0.499520 | 0.00048 | 0.30 |

![Coin tosses](plots/coin_tosses.png)

**Comment.**
- **Convergence.** fₙ converges to 0.5, as the **law of large numbers** says. But it
  converges **slowly and irregularly**:
  - after 100 tosses it is still at 0.62;
  - it touches 0.5 exactly at n = 1000 and then moves away again.
- **How fast.** The typical distance from 0.5 is the standard deviation of fₙ:
  σₙ = √(p(1−p)/n) = 1/(2√n). It shrinks only like **1/√n**: 100× more tosses make the error
  just 10× smaller. The deviations stay within about ±2σₙ (orange band), as expected (the
  2.4σ at n = 100 is a normal fluctuation).
- **What the law says.** The law of large numbers guarantees that the *fraction* converges. It
  does not say that the *number* of heads approaches n/2. At n = 100 000 there are
  48 fewer heads than half, and that difference typically grows like √n.

## 2. Monte Carlo estimate of π

Points (x, y) uniform in the unit square. The fraction falling inside the quarter circle
x² + y² < 1 estimates its area π/4, so π_N = 4·(points inside)/N.

**Two RNG streams.** x and y come from **two independent generators**, two Mersenne Twisters
with seeds 12345 and 12346, so the two coordinates of a point are not correlated.

**Error.** Each point is "inside" with probability p = π/4, so the theoretical error is
σ = 4√(p(1−p)/N) = √(π(4−π)/N) ≈ 1.64/√N. To measure the error, and not just one lucky or
unlucky run, every N is repeated R = 20 times, giving the RMS error:

| N | π (1st run) | Error | RMS error (20 runs) | σ theory |
|---|---|---|---|---|
| 10 | 3.2000000 | 5.8·10⁻² | 5.5·10⁻¹ | 5.2·10⁻¹ |
| 10² | 3.4400000 | 3.0·10⁻¹ | 1.6·10⁻¹ | 1.6·10⁻¹ |
| 10³ | 3.1280000 | −1.4·10⁻² | 4.4·10⁻² | 5.2·10⁻² |
| 10⁴ | 3.1388000 | −2.8·10⁻³ | 1.6·10⁻² | 1.6·10⁻² |
| 10⁵ | 3.1345600 | −7.0·10⁻³ | 4.8·10⁻³ | 5.2·10⁻³ |
| 10⁶ | 3.1422520 | 6.6·10⁻⁴ | 1.3·10⁻³ | 1.6·10⁻³ |
| 10⁷ | 3.1409968 | −6.0·10⁻⁴ | 5.5·10⁻⁴ | 5.2·10⁻⁴ |

**Fitted slope:** RMS error ∝ N^b with **b = −0.505** (theory −0.5).

![Monte Carlo pi](plots/monte_carlo_pi.png)

**Comment.**
- **The error decreases like 1/√N, independently of the dimension.** This is why Monte Carlo
  wins for high-dimensional integrals. It is also slow: every extra correct digit costs
  **100× more points**, and 10⁷ points give only ~3 correct decimals.
- **A single run is not a measure of the error.** At N = 10, the first run (3.2) is much closer
  to π than the typical error. At N = 10⁵, the error happens to be larger than at N = 10⁴. Only
  the RMS over many runs follows the 1/√N law.

## 3. Change of variables: Y = U²

For y ∈ (0, 1):
- the CDF is F_Y(y) = P(U² ≤ y) = P(U ≤ √y) = **√y**;
- so the density is f_Y(y) = **1/(2√y)**.

It diverges at y = 0 but is integrable. n = 10⁶ samples in 50 bins:

| | Sample | Exact |
|---|---|---|
| Mean E[Y] = 1/3 | 0.333114 | 0.333333 |
| Variance Var[Y] = 1/5 − 1/9 = 4/45 | 0.088850 | 0.088889 |
| χ² (49 degrees of freedom) | 41.0, p-value = 0.78 | |

![Y = U^2](plots/square_uniform.png)

**Comment.**
- **Agreement.** The histogram follows 1/(2√y). The log-log panel shows the power law
  y<sup>−1/2</sup> as a straight line.
- **The first bin.** Its height (7.09) is far from f_Y at the bin centre (f_Y(0.01) = 5.0),
  because the density changes a lot inside the bin. The correct comparison is with the
  *average* density over the bin, (√b − √a)/(b − a) = 7.07 for [0, 0.02) (dotted line).
- **The χ² test** compares observed and expected counts nᵢ = n(√bᵢ − √aᵢ) in all bins. A value
  close to the number of degrees of freedom (p = 0.78) means the histogram is consistent with
  f_Y.

## 4. Exponential distribution by inverse transform

**Inverse transform.** If F is a continuous CDF and U ~ Uniform(0,1), then Y = F⁻¹(U) has CDF
F, because P(F⁻¹(U) ≤ y) = P(U ≤ F(y)) = F(y). For the exponential, F(y) = 1 − e^(−λy), so

  **Y = −ln(1 − U)/λ**, with λ = 1.5.

Since U ∈ [0, 1), 1 − U ∈ (0, 1], so the logarithm is never ln 0. (−ln(U)/λ is equivalent,
because 1 − U is also uniform, but could hit U = 0.)

n = 10⁶ samples, 60 bins on [0, 6), plus an overflow bin y ≥ 6:

| | Sample | Exact |
|---|---|---|
| Mean 1/λ | 0.66618 | 0.66667 |
| Variance 1/λ² | 0.44389 | 0.44444 |
| Values ≥ 6 | 95 | 123.4 expected (n·e⁻⁹) |
| χ² (60 degrees of freedom) | 59.8, p-value = 0.48 | |

![Exponential](plots/exponential.png)

**Comment.** The histogram agrees with λe^(−λy) in both panels. On the log scale the
exponential is a straight line of slope −λ, and the agreement holds down to densities of 10⁻³.
Only in the tail, where each bin has a few hundred entries or fewer, do the statistical
fluctuations become visible: Poisson errors are ±√count.

## 5. Empirical CDF

The program uses the **same sample** as exercise 4 (same generator, seed and n). The empirical
CDF Fₙ(y) = (number of samples ≤ y)/n is computed by sorting the sample: it is a step function
jumping by 1/n at every value.

**Kolmogorov–Smirnov distance.** The distance between Fₙ and the exact
F(y) = 1 − e^(−λy) is measured by Dₙ = max_y |Fₙ(y) − F(y)|. If the sample really comes from
F, √n·Dₙ follows the Kolmogorov distribution, and **√n·Dₙ < 1.36 with 95% probability**.

| n (first n values) | Dₙ | √n·Dₙ |
|---|---|---|
| 10 | 0.195 | 0.62 |
| 100 | 0.159 | **1.59** |
| 1 000 | 0.030 | 0.95 |
| 10 000 | 0.0068 | 0.68 |
| 100 000 | 0.0019 | 0.61 |
| 1 000 000 | 0.00078 | 0.78 |

![Empirical CDF](plots/empirical_cdf.png)

**Comment.**
- **Left:** with 10⁶ points the empirical and exact CDFs are indistinguishable.
- **Centre:** the difference is at most 7.8·10⁻⁴.
- **Right:** Dₙ shrinks like 1/√n, the same √n law as in exercises 1 and 2. Most values of
  √n·Dₙ are well below 1.36.
- **The first 100 values exceed the 95% limit (1.59).** That is not an error. A correct sample
  exceeds the 95% limit in 5% of the cases, so 1 value out of 6 is not surprising, and these
  six checks are not independent anyway: each uses the first values of the same sample. With
  the full sample, √n·Dₙ = 0.78: consistent.
- **Why the empirical CDF is useful:** unlike a histogram, it needs no choice of bins.
