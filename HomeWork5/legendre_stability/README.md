# Stability of the Legendre recurrence

HomeWork05 asks to compute P<sub>ℓ</sub>(x), ℓ = 0 … 50, with the three-term recurrence
- forward (ℓ increasing),
- backward with Miller's algorithm (from L = 80, normalized with P<sub>0</sub> = 1),

and to compare both with high-precision reference values at x = 0.1, 0.5, 0.9, 0.99.

**Short answer: the homework's expectation is the wrong way round for P<sub>ℓ</sub>.** For
|x| < 1, the *forward* recurrence is accurate to ~10⁻¹⁶, and *Miller's* method gives wrong
results (errors of order 1). The professor's
["starting solution"](https://codimd.infn.it/s/ftKB8jdK1) reaches the same conclusion.

The story described in the homework (a dominant and a minimal solution, an unstable forward
recurrence, Miller as the fix) is real, but it happens for the **other** solution of the same
recurrence, **Q<sub>ℓ</sub>(x)**, at **x > 1**. The program shows both cases.

| File | Content |
|---|---|
| `legendre_stability.c` | Computes, compares, writes `output/legendre_<P\|Q>_x<value>.csv` |
| `plot_legendre.py` | Figures in `plots/` |
| `Makefile` | `make` (builds and prints the instructions), `make run`, `make plot`, `make help`, `make clean` |

## Build and run

```bash
make
./legendre_stability [l_max] [L] [x1 x2 ...]   # default: 50 80 0.1 0.5 0.9 0.99 1.1 1.5 3
make plot                                       # run + figures in plots/
```

For |x| < 1 the program studies P<sub>ℓ</sub>(x); for x > 1 it studies Q<sub>ℓ</sub>(x).

## 1. The three methods

All three use the recurrence (ℓ+1) f<sub>ℓ+1</sub> = (2ℓ+1) x f<sub>ℓ</sub> − ℓ f<sub>ℓ−1</sub>.
The starting values select which solution you get:

| Function | Starting values |
|---|---|
| P<sub>ℓ</sub> | P<sub>0</sub> = 1, P<sub>1</sub> = x |
| Q<sub>ℓ</sub> (x > 1) | Q<sub>0</sub> = ½ ln((x+1)/(x−1)), Q<sub>1</sub> = x Q<sub>0</sub> − 1 |

| Method | How |
|---|---|
| Forward | Recurrence upwards in `double` |
| Miller backward | f̃<sub>L+1</sub> = 0, f̃<sub>L</sub> = 1, recurrence downwards in `double`, then rescale so that f<sub>0</sub> = P<sub>0</sub> (or Q<sub>0</sub>) |
| **Reference** | `__float128` (113-bit mantissa, ~34 digits). For P<sub>ℓ</sub>: forward recurrence. For Q<sub>ℓ</sub>: Miller from ℓ<sub>max</sub> + 1000 |

**Checking the references.**
- The P<sub>ℓ</sub> reference agrees with **exact rational arithmetic** (Python `fractions`)
  to 10⁻²⁷, which is the number of digits written in the CSV.
- The Q<sub>ℓ</sub> reference agrees with a 120-digit computation to 10⁻³¹.
- The reference uses exactly the same (rounded) x as the double-precision code, so only the
  *algorithm's* error is measured.

**Why not GSL as the reference?** `gsl_sf_legendre_Pl_array` computes the same forward
recurrence in double precision, so its values are **bit-for-bit identical** to ours (column
"GSL =fw" in the program output). Measuring the error against GSL always gives exactly 0. That
was the problem with the first version of this exercise.

## 2–3. Errors (`make run`, l<sub>max</sub> = 50, L = 80)

Largest error over ℓ = 0 … 50:

| Function | Forward: abs | Forward: rel | Miller: abs | Miller: rel |
|---|---|---|---|---|
| P<sub>ℓ</sub>(0.1) | 2.8·10⁻¹⁷ | 6.9·10⁻¹⁵ | 2.5 | 2.8·10² |
| P<sub>ℓ</sub>(0.5) | 5.6·10⁻¹⁷ | 8.1·10⁻¹⁶ | 7.9 | 4.7·10¹ |
| P<sub>ℓ</sub>(0.9) | 2.5·10⁻¹⁶ | 5.6·10⁻¹⁵ | 0.12 | 1.6 |
| P<sub>ℓ</sub>(0.99) | 4.7·10⁻¹⁶ | 1.1·10⁻¹⁴ | 0.25 | 2.9 |
| Q<sub>ℓ</sub>(1.1) | 8.5·10⁻¹⁰ | **17** | 2.8·10⁻¹⁷ | 1.1·10⁻¹² |
| Q<sub>ℓ</sub>(1.5) | 6.3·10³ | **4.9·10²⁵** | 1.1·10⁻¹⁶ | 1.5·10⁻¹⁵ |
| Q<sub>ℓ</sub>(3) | 4.7·10²⁰ | **2.1·10⁶⁰** | 2.7·10⁻²⁰ | 4.6·10⁻¹⁶ |

| Claim in the homework | Result |
|---|---|
| Forward accurate far from 1 (x = 0.1) | $\color{green}{\textbf{O}}$ yes, ~10⁻¹⁶ |
| Forward unstable close to 1 (x = 0.99), errors of order 1 | $\color{red}{\textbf{X}}$ no: absolute error ≤ 5·10⁻¹⁶. Still accurate at ℓ = 1000 (`./legendre_stability 1000 1500 0.99`: 1.7·10⁻¹⁵) |
| Backward (Miller) accurate for all ℓ and x | $\color{red}{\textbf{X}}$ no: errors of order 1 for P<sub>ℓ</sub> at every x < 1 |
| *(same claims for Q<sub>ℓ</sub>, x > 1)* | $\color{green}{\textbf{O}}$ forward blows up, Miller accurate |

**Why are the relative errors of P<sub>ℓ</sub> sometimes 10⁻¹⁴ when the absolute ones are
10⁻¹⁶?** Near a zero of P<sub>ℓ</sub>(x), dividing by a tiny P<sub>ℓ</sub> inflates the
relative error. This is a property of the *function* (conditioning, section 6), not of the
algorithm. The dashed line κ<sub>ℓ</sub>(x)·ε in the plot is the error you would get just from
rounding x itself. The forward error follows it and does not grow with ℓ. It is above the line
only at a few points (at most ×20) where κ is accidentally almost 0, near the maxima of
P<sub>ℓ</sub>, where P′<sub>ℓ</sub> ≈ 0.

## 4. Connection to spherical harmonics

For m = 0,

  Y<sub>ℓ0</sub>(θ, φ) = √((2ℓ+1)/(4π)) · P<sub>ℓ</sub>(cos θ),

so Y<sub>ℓ0</sub> has exactly the same *relative* error as P<sub>ℓ</sub>(cos θ). The prefactor
is a known constant computed to full precision.

- Building Y<sub>ℓ0</sub> with the **forward** recurrence, as is usually done, is safe: relative
  error ~10⁻¹⁶ away from the zeros.
- Building it with **Miller's** recurrence would give completely wrong spherical harmonics:
  relative errors of 1–100 for ℓ ≥ 1.

## 5. Relative error vs ℓ

| Figure | Content |
|---|---|
| `plots/legendre_P_relative_error.png` | \|relative error\| of P<sub>ℓ</sub>(x), one panel per x: forward (blue) at the 10⁻¹⁶ level, Miller (red) at 10⁻¹–10² |
| `plots/legendre_Q_relative_error.png` | \|relative error\| of Q<sub>ℓ</sub>(x): forward grows as a straight line on the log scale (exponentially), Miller stays at 10⁻¹⁵ |
| `plots/legendre_values.png` | The values themselves: Miller's P<sub>ℓ</sub> visibly wrong, forward Q<sub>ℓ</sub> leaving the reference and growing |

![Relative error of P_l](plots/legendre_P_relative_error.png)

![Relative error of Q_l](plots/legendre_Q_relative_error.png)

![Values](plots/legendre_values.png)

## 6. Explanation

**Two solutions.** A three-term recurrence is a second-order difference equation, so it has
**two** independent solutions. For Bonnet's recurrence they are P<sub>ℓ</sub>(x) and
Q<sub>ℓ</sub>(x). In floating point, every step adds a rounding error of ~ε, which is a small
amount of *some* combination of both solutions. What happens to that error depends on how big
the two solutions are compared with each other:

| Region | P<sub>ℓ</sub> | Q<sub>ℓ</sub> | Which is dominant? |
|---|---|---|---|
| \|x\| < 1, x = cos θ | ≈ √(2/(πℓ sin θ)) cos((ℓ+½)θ − π/4) | similar amplitude, phase shifted | **neither**: both oscillate with amplitude ~ℓ<sup>−1/2</sup> |
| x > 1 | ~ λ<sup>ℓ</sup>/√ℓ, grows | ~ λ<sup>−ℓ</sup>/√ℓ, decays | P<sub>ℓ</sub> dominant, Q<sub>ℓ</sub> **minimal** |

where λ = x + √(x²−1) > 1.

- **Forward recurrence amplifies the dominant solution.** For x > 1, a rounding error of size
  ε·Q<sub>ℓ</sub> containing a bit of P<sub>ℓ</sub> grows like λ<sup>2ℓ</sup> *relative to*
  Q<sub>ℓ</sub>. The relative error of Q<sub>ℓ</sub> is therefore ≈ ε·λ<sup>2ℓ</sup>. This
  predicts 1.4·10²⁶ at x = 1.5 and 8·10⁶⁰ at x = 3, against the measured 4.9·10²⁵ and 2.1·10⁶⁰.
  For P<sub>ℓ</sub> at x > 1 there is no problem, because it is itself the dominant solution.
- **Backward recurrence (Miller's idea) suppresses it.** Going *down* in ℓ, the roles swap:
  Q<sub>ℓ</sub> grows and P<sub>ℓ</sub> shrinks. Whatever arbitrary start (f̃<sub>L</sub> = 1,
  f̃<sub>L+1</sub> = 0) is used, the P part dies out like λ<sup>−2(L−ℓ)</sup>, and after enough
  steps only Q<sub>ℓ</sub> remains, up to a constant fixed by the normalization at ℓ = 0.
  - At x = 1.1, with λ = 1.56 and L = 80, the leftover at ℓ = 50 is ≈ λ<sup>−60</sup> = 3·10⁻¹²
    (measured 1.1·10⁻¹²).
  - With L = 200 it disappears (10⁻¹⁵): `./legendre_stability 50 200 1.1`.
  - **The closer x is to 1, the larger L must be.**
- **Why Miller fails for P<sub>ℓ</sub> at |x| < 1.** Miller only works if the wanted solution is
  the *minimal* one, so that the backward recurrence suppresses everything else. For |x| < 1,
  P<sub>ℓ</sub> and Q<sub>ℓ</sub> have the same size, so nothing is suppressed in either
  direction. The backward recurrence returns an arbitrary combination
  a·P<sub>ℓ</sub> + b·Q<sub>ℓ</sub>, fixed by the arbitrary start. The single normalization
  f<sub>0</sub> = 1 fixes one constant, but cannot remove the Q<sub>ℓ</sub> part. The result
  therefore does not converge to P<sub>ℓ</sub>, however large L is. For the same reason, the
  forward recurrence has nothing to amplify, so it is **stable** for P<sub>ℓ</sub>: the error
  grows only slowly, roughly ∝ ℓ·ε.

**Conditioning vs stability.**
- **Conditioning** is a property of the *mathematical problem*: how much P<sub>ℓ</sub>(x)
  changes when x changes. The relative condition number is
  κ<sub>ℓ</sub>(x) = |x P′<sub>ℓ</sub>(x) / P<sub>ℓ</sub>(x)|. Here κ is moderate: the median
  over ℓ ≤ 50 is 2–142 for the four x values, and exactly at x = 1, κ = ℓ(ℓ+1)/2 ≈ 1275 at
  ℓ = 50. It is large only near the zeros of P<sub>ℓ</sub>, where no method can give a small
  *relative* error. So P<sub>ℓ</sub>(x) is well-conditioned near x = 1, as the homework says.
- **Stability** is a property of the *algorithm*: whether it adds more error than the
  conditioning makes unavoidable.
  - The forward recurrence for P<sub>ℓ</sub> stays at the κ·ε level and doesn't grow with ℓ:
    it is **stable**.
  - The forward recurrence for Q<sub>ℓ</sub> at x > 1 is **unstable**, with errors up to 10⁶⁰,
    even though Q<sub>ℓ</sub>(x) is a perfectly well-behaved, smoothly decaying function.
  - Miller's method is stable for Q<sub>ℓ</sub>, and simply *not applicable* to P<sub>ℓ</sub>.

  A well-conditioned problem can still be ruined by an unstable algorithm. Which direction is
  stable depends on which solution you want, not on the recurrence formula itself.
