# HomeWork 07: Random numbers

| Folder | Content |
|---|---|
| [`random_numbers/`](random_numbers/README.md) | Five C programs using the lecture notebook's LCG as the only random source, with their figures: coin tosses and the law of large numbers, Monte Carlo π with two streams, Y = U², exponential by inverse transform, empirical CDF, and a bonus on RANDU's planes |

**Main results:**
- The coin-toss fraction converges to 0.5 like 1/(2√n). The lowest bit of the LCG is useless:
  it just alternates 0, 1, 0, 1.
- The Monte Carlo error on π decreases like N^−0.500. The two "independent" streams are
  1.7·10⁹ steps apart on the same cycle of 2³² numbers.
- The histograms of U² and of the exponential pass a χ² test (p = 0.75 and 0.55).
- The empirical CDF of the exponential sample passes the Kolmogorov–Smirnov test
  (√n·Dₙ = 1.10 < 1.36).
- **Bonus:** all RANDU triples lie on 15 planes (9u<sub>n</sub> − 6u<sub>n+1</sub> + u<sub>n+2</sub>
  is always an integer).
