# HomeWork05: Stability of Legendre polynomials

## The assignment

For order m = 0 and degrees up to ℓ<sub>max</sub> = 50, compare numerical ways of computing
the Legendre polynomials P<sub>ℓ</sub>(x) (which are Y<sub>ℓ0</sub> up to normalization) at
**x = 0.1, 0.5, 0.9, 0.99**, and check whether the expected stability claim is actually true.
The three methods asked for:

1. **Direct three-term forward recurrence** (Bonnet):
   P<sub>0</sub> = 1, P<sub>1</sub> = x,
   P<sub>ℓ+1</sub>(x) = (2ℓ+1)/(ℓ+1) · x P<sub>ℓ</sub>(x) − ℓ/(ℓ+1) · P<sub>ℓ−1</sub>(x),  ℓ ≥ 1.
2. **High-precision reference values**, treated as exact when measuring the double-precision
   errors.
3. **Backward recurrence experiment** (Miller style): solve the recurrence for
   P<sub>ℓ−1</sub>,
   P<sub>ℓ−1</sub>(x) = (2ℓ+1)/ℓ · x P<sub>ℓ</sub>(x) − (ℓ+1)/ℓ · P<sub>ℓ+1</sub>(x),  ℓ = L, L−1, …, 1,
   start from P̃<sub>L+1</sub> = 0, P̃<sub>L</sub> = 1, propagate backward, rescale with
   c = 1/P̃<sub>0</sub>, P<sub>ℓ</sub><sup>(back)</sup> = c·P̃<sub>ℓ</sub>, and check whether this
   recovers P<sub>ℓ</sub>(x).

List of actions: (1) implement the methods; (2) for ℓ = 0 … 50 compute absolute and relative
errors; (3) plot the relative error versus ℓ for each x, and inspect the absolute error;
(4) explain the numerical behaviour; (5) connect the result to the spherical harmonics.

The sections below follow the list of actions.

| File | Content |
|---|---|
| `legendre_stability.c` | Computes the three methods and the errors, writes `output/legendre_P_x<value>.csv` |
| `plot_legendre.py` | Figures in `plots/` |
| `Makefile` | `make` (builds and prints the instructions), `make run`, `make plot`, `make help`, `make clean` |

## Build and run

```bash
make
./legendre_stability [l_max] [L] [x1 x2 ...]   # default: 50 80 0.1 0.5 0.9 0.99
make plot                                       # run + figures in plots/
```

## 1. The methods

| Method | How it is computed |
|---|---|
| **1. Forward** | Bonnet's recurrence upwards from P<sub>0</sub> = 1, P<sub>1</sub> = x, in `double` |
| **2. Reference** | The same forward recurrence in quad precision, `__float128` (113-bit mantissa, about 34 decimal digits) |
| **3. Backward (Miller style)** | P̃<sub>L+1</sub> = 0, P̃<sub>L</sub> = 1 with L = 80, the backward formula down to ℓ = 0 in `double`, then P<sub>ℓ</sub><sup>(back)</sup> = P̃<sub>ℓ</sub>/P̃<sub>0</sub> |

**The reference is reliable.**
- A quad-precision forward recurrence is "a high-precision implementation of the recurrence",
  as allowed by the assignment. With ~34 digits it is about 10¹⁸ times more precise than a
  `double`, so it can be treated as exact.
- It agrees with **exact rational arithmetic** (Python `fractions`) to 10⁻²⁷, the number of
  digits written in the CSV files.
- It uses exactly the same (rounded) x as the double-precision code, so only the error of the
  *algorithm* is measured.
- GSL (`gsl_sf_legendre_Pl_array`) could **not** be the reference: it computes the same forward
  recurrence in `double`, so its values are bit-for-bit identical to ours (column "GSL =fw" in
  the program output) and the measured error would always be exactly 0.

## 2. Absolute and relative errors (ℓ = 0 … 50)

`make run` (ℓ<sub>max</sub> = 50, L = 80). Largest error over ℓ = 0 … 50:

| x | Forward: abs | Forward: rel | Backward: abs | Backward: rel |
|---|---|---|---|---|
| 0.1 | 2.8·10⁻¹⁷ | 6.9·10⁻¹⁵ | 2.5 | 2.8·10² |
| 0.5 | 5.6·10⁻¹⁷ | 8.1·10⁻¹⁶ | 7.9 | 4.7·10¹ |
| 0.9 | 2.5·10⁻¹⁶ | 5.6·10⁻¹⁵ | 0.12 | 1.6 |
| 0.99 | 4.7·10⁻¹⁶ | 1.1·10⁻¹⁴ | 0.25 | 2.9 |

- **Forward:** absolute errors ≤ 5·10⁻¹⁶ at every x, also close to 1 (x = 0.99). It stays
  accurate for much larger ℓ too: `./legendre_stability 1000 1500 0.99` gives 1.7·10⁻¹⁵ at
  ℓ = 1000.
- **Backward:** errors of order 1 (relative errors from 1.6 to 280): the result is simply not
  P<sub>ℓ</sub>(x). A larger starting point does not help: with L = 200 or L = 1000 the errors
  at x = 0.5 and 0.99 are still between 0.2 and 15 (e.g. `./legendre_stability 50 1000 0.5`).
- **Why are the relative errors sometimes 10⁻¹⁴ when the absolute ones are 10⁻¹⁶?** Near a zero
  of P<sub>ℓ</sub>(x), dividing by a tiny P<sub>ℓ</sub> inflates the relative error. This is a
  property of the *function* (its conditioning, section 4), not of the algorithm. That's why the
  absolute error must be inspected too.

## 3. Plots

`make plot` writes three figures in `plots/`:

![Relative error of P_l](plots/legendre_P_relative_error.png)

*`plots/legendre_P_relative_error.png`: |relative error| of P<sub>ℓ</sub>(x) versus ℓ, one
panel per x. Forward (blue) at the 10⁻¹⁶ level, backward (red) at 10⁻¹–10². The dashed grey
line κ<sub>ℓ</sub>(x)·ε is the error caused just by rounding x itself (conditioning,
section 4).*

![Absolute error of P_l](plots/legendre_P_absolute_error.png)

*`plots/legendre_P_absolute_error.png`: |absolute error| of P<sub>ℓ</sub>(x) versus ℓ. The
forward error never exceeds a few times ε; the backward error is of order 0.1–10 for every
ℓ ≥ 1.*

![Values](plots/legendre_values.png)

*`plots/legendre_values.png`: the values themselves. The forward values lie on the reference;
the backward ones are a different oscillating sequence that only agrees with P<sub>ℓ</sub> at
ℓ = 0, where it is forced to by the rescaling.*

## 4. Explanation of the numerical behaviour

**A three-term recurrence has two independent solutions.** Bonnet's recurrence is a
second-order difference equation, so its solutions form a two-dimensional space: every solution
is a·P<sub>ℓ</sub>(x) + b·Q<sub>ℓ</sub>(x), where Q<sub>ℓ</sub> is the Legendre function of
the second kind. The two starting values choose a and b. In floating point every step adds a
rounding error of ~ε, which is a small amount of *some* combination of both solutions. What
happens to it depends on how big the two solutions are compared with each other. For
|x| < 1, with x = cos θ:

P<sub>ℓ</sub>(cos θ) ≈ √(2/(πℓ sin θ)) · cos((ℓ+½)θ − π/4),

and Q<sub>ℓ</sub> has the same form with a shifted phase: **both oscillate with the same
amplitude ~ℓ<sup>−1/2</sup>**. Neither solution dominates the other.

**Miller's algorithm is useful when the wanted solution is minimal in the direction of
backward propagation.** If one solution is *minimal* (it decreases going up in ℓ while the
other grows), then going *down* in ℓ the wanted solution is the one that grows. Whatever the
arbitrary start (P̃<sub>L</sub> = 1, P̃<sub>L+1</sub> = 0), the backward recurrence quickly
forgets the unwanted part, and after enough steps only the wanted solution is left, up to a
constant that the rescaling c = 1/P̃<sub>0</sub> fixes.

**For P<sub>ℓ</sub>(x) on |x| < 1 there is no dominant/minimal separation.** P<sub>ℓ</sub> and
Q<sub>ℓ</sub> have the same size, so nothing is suppressed in either direction. The backward
recurrence returns an arbitrary combination a·P<sub>ℓ</sub> + b·Q<sub>ℓ</sub>, fixed by the
arbitrary start. The rescaling with c = 1/P̃<sub>0</sub> fixes only one constant (the value at
ℓ = 0) and cannot remove the Q<sub>ℓ</sub> part, so the result is not P<sub>ℓ</sub>, however
large L is (section 2).

**Therefore the forward recurrence is stable.** For the same reason, going forward there is no
dominant solution to amplify the rounding errors: they stay at the ε level and grow only
slowly, roughly ∝ ℓ·ε. This is what the errors and the plots show.

**Conditioning vs stability.**
- **Conditioning** is a property of the *mathematical problem*: how much P<sub>ℓ</sub>(x)
  changes when x changes. The relative condition number is
  κ<sub>ℓ</sub>(x) = |x P′<sub>ℓ</sub>(x) / P<sub>ℓ</sub>(x)|. Here κ is moderate (median over
  ℓ ≤ 50 between 2 and 142 for the four x values; exactly at x = 1, κ = ℓ(ℓ+1)/2 ≈ 1275 at
  ℓ = 50). It is large only near the zeros of P<sub>ℓ</sub>, where no method can give a small
  *relative* error.
- **Stability** is a property of the *algorithm*: whether it adds more error than the
  conditioning makes unavoidable. The forward recurrence stays at the κ·ε level (the dashed line
  in the relative-error plot; the forward error is above it only at a few points, at most ×20,
  where κ is accidentally almost 0 near the maxima of P<sub>ℓ</sub>) and does not grow with ℓ:
  it is **stable**. The backward procedure is not an unstable version of the same computation:
  it computes a *different* solution of the recurrence.

## 5. Connection to spherical harmonics

For m = 0,

  Y<sub>ℓ0</sub>(θ, φ) = √((2ℓ+1)/(4π)) · P<sub>ℓ</sub>(cos θ),

so an error in P<sub>ℓ</sub>(cos θ) goes directly into Y<sub>ℓ0</sub>: the prefactor is a
known constant computed to full precision, and Y<sub>ℓ0</sub> has exactly the same *relative*
error as P<sub>ℓ</sub>(cos θ).
- Building Y<sub>ℓ0</sub> with the **forward** recurrence, as is usually done, is safe:
  relative error ~10⁻¹⁶ away from the zeros.
- Building it with the **Miller-style backward** procedure would give completely wrong
  spherical harmonics: relative errors of 1–100 for ℓ ≥ 1.

## Conclusion

For x = 0.1, 0.5, 0.9, 0.99 and ℓ<sub>max</sub> = 50:
- the ordinary Legendre **forward** recurrence is **stable** for all the tested cases, with
  absolute errors ≤ 5·10⁻¹⁶ against the high-precision reference;
- the arbitrary Miller-style **backward** procedure does **not** compute P<sub>ℓ</sub>(x) for
  |x| < 1: its errors are of order 1, whatever the starting point L;
- **Miller's algorithm is useful only where the wanted solution is genuinely minimal**, which
  is not the case for P<sub>ℓ</sub>(x) on |x| < 1, where P<sub>ℓ</sub> and Q<sub>ℓ</sub> have
  the same size.
