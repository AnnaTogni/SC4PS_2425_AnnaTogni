# HomeWork 05: Pointers and floating-point numbers

| Folder | Exercise | Main result |
|---|---|---|
| [`floating_point_sum/`](floating_point_sum/README.md) | Lecture exercise: sum [1, 10¹⁶, −10¹⁶, −0.5] with a for loop, GSL and Kahan | All three give −0.5 instead of 0.5. Kahan fails because the new term is larger than the running sum; Neumaier's variant and an exact summation give 0.5 |
| [`legendre_stability/`](legendre_stability/README.md) | HomeWork05: stability of Legendre polynomials (forward vs Miller-style backward recurrence vs high-precision reference) | For P<sub>ℓ</sub>(x), \|x\| < 1, the forward recurrence is stable (errors ~10⁻¹⁶); the Miller-style backward recurrence does not recover P<sub>ℓ</sub> (errors of order 1), because there is no minimal solution to select. Additional tests: Miller works for Q<sub>ℓ</sub>(x), x > 1, where the forward recurrence is unstable; a backward recurrence with two conditions (P<sub>0</sub> and P<sub>1</sub>) recovers P<sub>ℓ</sub> |

Each folder has a `Makefile`: `make` builds the program and prints how to use it, and
`make help` prints the instructions again.

Requirements: `gcc`, GNU GSL, `libquadmath` (part of gcc), Python 3 with `matplotlib` for the
Legendre plots.
