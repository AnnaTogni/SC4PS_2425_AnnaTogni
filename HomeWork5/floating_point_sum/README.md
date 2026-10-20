# Floating-point summation

Exercise from the floating-point lecture: sum the vector **[1.0, 1.0e16, −1.0e16, −0.5]**
(exact sum **0.5**)
- with a for loop,
- with a library function (GSL),
- with Kahan's algorithm,

and explain why the results differ.

`floating_point_sum.c` adds five more methods to understand what is going on.

| File | Content |
|---|---|
| `floating_point_sum.c` | The program |
| `Makefile` | `make` (builds and prints the instructions), `make run`, `make help`, `make clean` |

## Build and run

```bash
make
./floating_point_sum                  # the lecture's vector
./floating_point_sum 1e16 1 -1e16     # any list of numbers
```

## Results for [1, 1e16, −1e16, −0.5]

| Method | Result | Correct? |
|---|---|---|
| 1. for loop | −0.5 | $\color{red}{\textbf{X}}$ |
| 2. GSL `gsl_vector_sum` | −0.5 | $\color{red}{\textbf{X}}$ |
| 3. Kahan | −0.5 | $\color{red}{\textbf{X}}$ |
| 4. Neumaier (improved Kahan) | **0.5** | $\color{green}{\textbf{O}}$ |
| 5. Pairwise | 0 | $\color{red}{\textbf{X}}$ |
| 6. Sorted by \|x\|, then loop | 0 | $\color{red}{\textbf{X}}$ |
| 7. Loop in `long double` | **0.5** | $\color{green}{\textbf{O}}$ |
| 8. Exact sum (Shewchuk, as Python's `math.fsum`) | **0.5** | $\color{green}{\textbf{O}}$ |

## Why are the results different?

**The 1 disappears.** Near 10¹⁶, doubles are spaced by 2: a double has 53 bits of mantissa,
and 2⁵³ ≈ 9·10¹⁵. So 10¹⁶ + 1 is exactly halfway between two doubles, and it rounds (to even)
to 10¹⁶. Floating-point addition is **not associative**: (1 + 10¹⁶) − 10¹⁶ = 0, but
1 + (10¹⁶ − 10¹⁶) = 1.

- **For loop and GSL:** 1 + 10¹⁶ → 10¹⁶, then − 10¹⁶ → 0, then − 0.5 → **−0.5**.
  `gsl_vector_sum` is a plain loop, so it gives the same result.
- **Kahan:** the compensation c = (t − s) − y recovers the lost low-order part *only if
  |s| ≥ |y|* (the running sum is larger than the new term). At step 2, the new term (10¹⁶) is
  much larger than the sum (1), so t − s = 10¹⁶ − 1 is itself rounded to 10¹⁶, and c becomes 0.
  The lost 1 is never recorded:

  | Step | Add | Sum s | Compensation c |
  |---|---|---|---|
  | 1 | 1 | 1 | 0 |
  | 2 | 1e16 | 1e16 | **0** (the 1 is lost) |
  | 3 | −1e16 | 0 | 0 |
  | 4 | −0.5 | −0.5 | 0 |
- **Neumaier:** computes the lost part from whichever operand is *larger*, so it works in both
  cases. At step 2 it records c = 1 and adds it back at the end: −0.5 + 1 = **0.5**.
- **Pairwise:** (1 + 10¹⁶) + (−10¹⁶ − 0.5) = 10¹⁶ − 10¹⁶ = **0**. Both the 1 and the 0.5 are
  lost, because each one is added to a number of size 10¹⁶.
- **Sorted by \|x\|:** the order is −0.5, 1, 10¹⁶, −10¹⁶. The partial sum 0.5 is added to 10¹⁶
  and lost, so the result is **0**. Sorting helps for many numbers of similar sign, but not when
  large terms cancel.
- **long double:** 64 bits of mantissa. Near 10¹⁶ the spacing is 1/1024, so 10¹⁶ + 1 is exact and
  the result is **0.5**. More precision only moves the problem further away: try
  `./floating_point_sum 1 1e20 -1e20 -0.5`.
- **Exact (Shewchuk):** keeps the running sum *exactly*, as a list of non-overlapping doubles,
  using error-free transformations (a + b = hi + lo exactly). It always gives the correctly
  rounded result.

**Takeaways.**
- The result of a floating-point sum depends on the **order** of the operations.
- Compensated algorithms can fix it, but **Kahan is not enough** when a term is larger than the
  running sum. Neumaier, or an exact summation, is needed.
- One more surprise: `./floating_point_sum 0.1 0.2 0.3 -0.6` gives **2.8·10⁻¹⁷**, not 0, even
  with the exact method. The doubles closest to 0.1, 0.2, 0.3 and 0.6 are not exactly those
  decimal numbers, and their exact sum really is 2.8·10⁻¹⁷.
