# HomeWork05: Stability of Legendre Polynomial

This solution implements the three requested computations for `P_l(x)` with `m = 0`:

- forward three-term recurrence in `double`
- backward Miller recurrence in `double`
- a quad-precision reference recurrence for error measurements
- a GSL reference as an additional cross-check

## Build

```bash
make
```

## Run

```bash
make run
```

## Plot

```bash
make plot
```

This generates PNG comparison plots using Python/Matplotlib in the `plots_python/` folder.

Aliases also work:

```bash
make python_plot
make python
make pplot
make python pplot
```

The program generates one CSV file per test value:

- `legendre_x_0p1.csv`
- `legendre_x_0p5.csv`
- `legendre_x_0p9.csv`
- `legendre_x_0p99.csv`

Each CSV includes:

- `l`
- the forward result
- the Miller backward result
- the quad-precision reference
- the GSL reference
- signed and relative errors for both recurrences

The CSV columns are, in order:
1. `l`
2. `x`
3. `forward_double`
4. `backward_miller`
5. `reference_quad`
6. `reference_gsl`
7. `err_forward`
8. `rel_err_forward`
9. `err_backward`
10. `rel_err_backward`


These files are ready to plot with Python or inspect in a spreadsheet.

## Exercise Points 3-6

### 3. Loss of significance and instability

The Python plots now include:

- the values of `P_l(x)` computed by the three methods
- the signed error
- the signed relative error
- the relative error magnitude `|rel_err|` on a logarithmic scale

This last panel is the one to use for the exercise statement “plot relative error vs. `l`”.
The full figure now contains four panels:

- `P_l(x)` values
- signed error
- signed relative error
- absolute relative error

With the current implementation, the forward recurrence is indeed reasonably accurate for values such as `x = 0.1`, and the quad-precision reference lets you see small but nonzero deviations.

However, for plain Legendre polynomials `P_l(x)` with `m = 0`, the forward recurrence does **not** blow up dramatically on the default range `l <= 50`. To expose stronger finite-precision effects you should run more demanding cases, for example:

```bash
./legendre_stability 1000 1400 0.1 0.99 0.9999 0.99999
python3 plot_legendre.py
```

This makes the relative-error panel much more informative near `x = 1`.

Important note:

- the present Miller implementation is the direct textbook backward recursion with simple normalization at `P_0 = 1`
- for this `m = 0` test it does **not** yet outperform the forward recurrence uniformly
- so the code now documents the numerical behavior honestly instead of forcing the expected conclusion

### 4. Connection to spherical harmonics

For `m = 0`,

$$
Y_{\ell 0}(\theta,\phi) = \sqrt{\frac{2\ell+1}{4\pi}}\,P_\ell(\cos\theta).
$$

So any numerical error in `P_l(x)` is inherited directly by `Y_{l0}` up to the prefactor

$$
\sqrt{\frac{2\ell+1}{4\pi}}.
$$

That means an unstable recurrence for `P_l` contaminates the corresponding spherical harmonics as well.

### 5. Relative error vs. `l`

Use the fourth panel of each Python plot:

- it shows `|rel_err_forward|` and `|rel_err_backward|`
- it is the clearest panel for checking whether either method grows with `l`
- this is the panel that corresponds most closely to the exercise request

### 6. Short numerical explanation

Three-term recurrences usually admit two independent solutions:

- a dominant solution, which grows faster under the recurrence
- a minimal solution, which is easily polluted by rounding errors

Forward recurrence is vulnerable because tiny roundoff components in the dominant solution can be amplified step after step. Even if the exact target solution is well behaved, the algorithm can drift toward the dominant mode.

Backward recurrence (Miller’s idea) reverses the direction. In favorable problems, the unwanted dominant component is damped while the desired minimal solution is recovered after a final normalization.

This is the key distinction between conditioning and stability:

- the mathematical quantity `P_l(x)` can still be well conditioned
- but a particular finite-precision algorithm for computing it can be unstable

So the issue is often not “the function is impossible to evaluate”, but rather “this recurrence direction is a bad numerical algorithm for that regime”.

## Optional arguments

You can override the default setup:

```bash
./legendre_stability [l_max] [miller_L] [x1 x2 ...]
```

Example:

```bash
./legendre_stability 200 400 0.9 0.99 0.999
```

To do the plots using the non default set up:

```bash
mkdir -p plots_python
MPLCONFIGDIR=/tmp/mplconfig python3 plot_legendre.py
```
