# Floating-point summation

Exercise from the floating-point lecture: define the vector
$$v = [1.0,\ 10^{16},\ -10^{16},\ -0.5],$$

whose exact sum is **0.5**, and compute its sum in three
different ways:
1. with a `for` loop;
2. with the sum function of a library (GSL, `gsl_vector_sum`);
3. with Kahan's summation algorithm;

then compare the results.

| File | Content |
|---|---|
| `floating_point_sum.c` | The program: the 3 required methods, the 5 additional ones, and a step-by-step trace of Kahan and Neumaier |
| `Makefile` | `make` (builds and prints the instructions), `make run`, `make help`, `make clean` |

## Build and run

```bash
make
./floating_point_sum                  # the lecture's vector
./floating_point_sum 1e16 1 -1e16     # any list of numbers
```

## Results for [1, 1e16, −1e16, −0.5]

| Method | | Result | |
|---|---|---|---|
| 1. `for` loop | required | −0.5 | $\color{red}{\textbf{wrong}}$ |
| 2. GSL `gsl_vector_sum` | required | −0.5 | $\color{red}{\textbf{wrong}}$ |
| 3. Kahan | required | −0.5 | $\color{red}{\textbf{wrong}}$ |
| 4. Neumaier (improved Kahan) | additional | **0.5** | $\color{green}{\textbf{correct}}$ |
| 5. Pairwise | additional | 0 | $\color{red}{\textbf{wrong}}$ |
| 6. Sorted by \|x\|, then loop | additional | 0 | $\color{red}{\textbf{wrong}}$ |
| 7. Loop in `long double` | additional | **0.5** | $\color{green}{\textbf{correct}}$ |
| 8. Exact sum (Shewchuk, as Python's `math.fsum`) | additional | **0.5** | $\color{green}{\textbf{correct}}$ |

## The key fact: the 1 disappears

A `double` has 53 bits of mantissa (52 stored + 1 implicit). Between $2^{53} \approx 9.0 \cdot 10^{15}$
and $2^{54} \approx 1.8 \cdot 10^{16}$ the doubles are therefore **2 apart**: $10^{16}$ and
$10^{16} + 2$ are doubles, $10^{16} + 1$ is not. $10^{16} + 1$ lies exactly halfway, and IEEE 754
rounds it (to the even neighbour) back to $10^{16}$. The 1 is lost.

As a consequence, floating-point addition is **not associative**: the result depends on the
order of the operations:

$$(1 + 10^{16}) - 10^{16} = 10^{16} - 10^{16} = 0, \qquad \text{but} \qquad 1 + (10^{16} - 10^{16}) = 1 + 0 = 1.$$

## The three required methods

### 1. `for` loop

```c
double s = 0.0;
for (size_t i = 0; i < n; i++) s += v[i];
```

The elements are added left to right: $1 + 10^{16} \to 10^{16}$ (the 1 is lost), then
$-10^{16} \to 0$, then $-0.5 \to -0.5$.

### 2. GSL `gsl_vector_sum`

The numbers are wrapped in a GSL vector (`gsl_vector_view_array`, no copy) and summed with
`gsl_vector_sum`. Internally this function is a plain loop over the elements, so it makes the
same roundings in the same order: **−0.5**. A library function is not automatically more
accurate.

### 3. Kahan compensated summation

**Idea** (Kahan, 1965): keep, next to the running sum `s`, a variable `c` with the small part
lost by the rounding of the last addition, and put it back into the next term.

```c
double s = 0.0, c = 0.0;
for (size_t i = 0; i < n; i++) {
    double y = v[i] - c;   // the new term, corrected by what was lost before
    double t = s + y;      // the new sum (rounded)
    c = (t - s) - y;       // (what was really added) - (what we wanted to add) = the lost part
    s = t;
}
```

For many numbers of similar size this makes the error independent of $n$. Here it does not help:
the formula `c = (t − s) − y` is exact **only if** $\lvert s \rvert \ge \lvert y \rvert$, i.e. if the running sum is larger
than the new term. At step 2 the new term ($10^{16}$) is much larger than the sum (1): `t − s` =
$10^{16} - 1$ is itself rounded to $10^{16}$, so `c` = 0 and the lost 1 is never recorded. The program
prints this step by step:

| Step | Add | Sum s | Compensation c |
|---|---|---|---|
| 1 | 1 | 1 | 0 |
| 2 | 1e16 | 1e16 | **0** (the 1 is lost) |
| 3 | −1e16 | 0 | 0 |
| 4 | −0.5 | −0.5 | 0 |

Result: **−0.5**.

## Five additional methods

These are not part of the exercise. Each one tests one idea about *why* the required methods
fail and *how* the sum can be computed correctly.

### 4. Neumaier: does a better compensation fix Kahan? Yes

**Idea** (Neumaier, 1974): compute the lost part from whichever of the two operands is
**larger**, so that the formula is exact in both cases. The corrections are accumulated in `c`
and added only once, at the end.

```c
double t = s + v[i];
if (fabs(s) >= fabs(v[i])) c += (s - t) + v[i];   // the sum is larger: as in Kahan
else                       c += (v[i] - t) + s;   // the new term is larger
s = t;
...
return s + c;
```

At step 2, $\lvert s \rvert = 1 < \lvert v_i \rvert = 10^{16}$, so `c += (1e16 − 1e16) + 1 = 1`:
the lost 1 is recorded. At the end: $-0.5 + 1 = 0.5$. Correct.

| Step | Add | Sum s | Correction c |
|---|---|---|---|
| 1 | 1 | 1 | 0 |
| 2 | 1e16 | 1e16 | **1** (the 1 is saved) |
| 3 | −1e16 | 0 | 1 |
| 4 | −0.5 | −0.5 | 1 |

### 5. Pairwise summation: does a different order help? Not here

**Idea:** split the vector in two halves, sum each half recursively, then add the two results.
The rounding error grows like $\log_2 n$ instead of $n$, which is why it is used, for example, by
NumPy's `sum`.

Here:

$$(1 + 10^{16}) + (-10^{16} - 0.5) = 10^{16} + (-10^{16}) = 0.$$

Both the 1 and the $-0.5$ are added to a number of size $10^{16}$ and are lost.

### 6. Sort by $\lvert x \rvert$, then loop: does "small numbers first" help? Not here

**Idea:** a classic rule is to add the numbers from the smallest to the largest in absolute
value, so that the small ones are accumulated before they meet the big ones.

The order becomes $-0.5,\ 1,\ 10^{16},\ -10^{16}$: $-0.5 + 1 = 0.5$, then
$0.5 + 10^{16} \to 10^{16}$ (the 0.5 is lost), then $-10^{16} \to 0$. The rule works when all the numbers have the same sign, but not when
large terms cancel each other.

### 7. Loop in `long double`: does more precision help? Yes, but only here

**Idea:** the same plain loop, with the running sum in `long double`. On x86 this is the 80-bit
extended format, with 64 bits of mantissa instead of 53.

Near $10^{16}$ the spacing becomes $2^{53-63} = 1/1024$, so $10^{16} + 1$ is stored exactly and the result is
**0.5**. More precision only moves the problem further away: with larger numbers it comes back,
e.g. `./floating_point_sum 1 1e20 -1e20 -0.5` fails again.

### 8. Exact summation (Shewchuk): what is the true answer? 0.5

**Idea** (Shewchuk, 1997; the algorithm behind Python's `math.fsum`): keep the running sum
**exactly**, as a short list of doubles that do not overlap ("partials"). Each new number is
added to the list with an *error-free transformation*: for two doubles $a$ and $b$ with $\lvert a \rvert \ge \lvert b \rvert$,

$$\text{hi} = \mathrm{fl}(a + b), \qquad \text{lo} = b - (\text{hi} - a)$$

gives two doubles with $\text{hi} + \text{lo} = a + b$ **exactly** (hi is the rounded sum, lo the
rounding error). Only at the end are the partials added into one double, correctly rounded.

Result: **0.5**, the correctly rounded value of the exact sum. It is used here as the reference:
it confirms that 0.5 is the right answer, and that the other methods are wrong because of
rounding, not because of a bug.

## Summary

- The result of a floating-point sum depends on the **order** of the operations, because each
  addition is rounded and addition is not associative.
- The three required methods all lose the 1 when it is added to $10^{16}$, so they give $-0.5$.
- **Kahan's compensation is not enough** when a term is larger than the running sum; Neumaier's
  version fixes exactly this case.
- Changing the order (pairwise, sorting) does not help when large terms cancel.
- More precision (`long double`) helps only until the numbers become larger.
- An exact summation (Shewchuk / `math.fsum`) always gives the correctly rounded result.

## References

- The lecture note on floating-point numbers: <https://codimd.infn.it/s/QGAx8wKwd>
- Kahan summation algorithm (Wikipedia, linked in the exercise):
  <https://en.wikipedia.org/wiki/Kahan_summation_algorithm>
- W. Kahan, *Further remarks on reducing truncation errors*, Communications of the ACM **8**(1),
  40 (1965). The original Kahan algorithm.
- A. Neumaier, *Rundungsfehleranalyse einiger Verfahren zur Summation endlicher Summen*,
  ZAMM (Zeitschrift für Angewandte Mathematik und Mechanik) **54**(1), 39–51 (1974). The
  improved version (method 4).
- J. R. Shewchuk, *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric
  Predicates*, Discrete & Computational Geometry **18**(3), 305–363 (1997). Exact sums with
  error-free transformations (method 8).
- T. J. Dekker, *A floating-point technique for extending the available precision*, Numerische
  Mathematik **18**, 224–242 (1971). The error-free transformation hi + lo = a + b.
- N. J. Higham, *The accuracy of floating point summation*, SIAM Journal on Scientific Computing
  **14**(4), 783–799 (1993). Comparison of loop, pairwise, sorted and compensated summation.
- N. J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM (2002),
  chapter 4 "Summation".
- D. Goldberg, *What every computer scientist should know about floating-point arithmetic*,
  ACM Computing Surveys **23**(1), 5–48 (1991).
- IEEE Standard for Floating-Point Arithmetic, IEEE 754-2019.
- GSL (GNU Scientific Library) reference manual, *Vectors* (`gsl_vector_sum`):
  <https://www.gnu.org/software/gsl/doc/html/vectors.html>
- Python documentation, `math.fsum`: <https://docs.python.org/3/library/math.html#math.fsum>
