# HomeWork 07: Random numbers

| Folder | Content |
|---|---|
| [`random_numbers/`](random_numbers/README.md) | Five C programs (GSL Mersenne Twister) and their figures: coin tosses and the law of large numbers, Monte Carlo π, Y = U², exponential by inverse transform, empirical CDF |

**Main results:**
- The coin-toss fraction converges to 0.5 like 1/(2√n).
- The Monte Carlo error on π decreases like N^−0.505 (theory N^−½).
- The histograms of U² and of the exponential pass a χ² test (p = 0.78 and 0.48).
- The empirical CDF of 10⁶ exponential samples differs from the exact CDF by at most 7.8·10⁻⁴
  (Kolmogorov–Smirnov √n·Dₙ = 0.78 < 1.36).
