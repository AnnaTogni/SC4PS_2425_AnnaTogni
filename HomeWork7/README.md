# HomeWork 07: Random numbers

| Folder | Content |
|---|---|
| [`random_numbers/`](random_numbers/README.md) | Five C programs, one per exercise, using the lecture notebook's LCG as the only random source, with their figures: coin tosses and the law of large numbers, Monte Carlo π with two streams, Y = U², exponential by inverse transform, empirical CDF |

**Main results:**
- The fraction of heads converges to 0.5, slowly: its distance from 0.5 shrinks like 1/(2√n).
- The Monte Carlo error on π decreases like N^−0.500.
- The histograms of U² and of the exponential sample agree with the analytic densities
  1/(2√y) and λe^(−λy).
- The empirical CDF of the exponential sample agrees with the exact CDF 1 − e^(−λy) to within
  2.5·10⁻³.
