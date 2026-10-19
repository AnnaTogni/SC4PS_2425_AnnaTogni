# HomeWork 02: vector sum d = a·x + y

Numerical computing's "Hello World": `NumericalSum_v2.c` computes **d = a·x + y**.
- **a** is a scalar.
- **x** and **y** are vectors of dimension N in which all elements are equal.

The program then checks that every element of **d** has the expected value.

| File | Content |
|---|---|
| `NumericalSum_v2.c` | Solution of the homework (current version) |
| `Makefile` | Builds `numericalsum_v2` |
| `NumericalSum.c` | First version, kept for reference |
| `CHANGES.md` | What changed between the two versions, and why |

## Build and run

```bash
make
./numericalsum_v2 N a x y [expected]
```

| Argument | Meaning |
|---|---|
| `N` | Vector dimension. Scientific notation is accepted (`1e6`, `1e8`). |
| `a` | Scalar |
| `x`, `y` | Value of every element of x and of y |
| `expected` | *(optional)* Value every `d[i]` should be equal to |

Example (point 4 of the assignment):
```bash
./numericalsum_v2 1e8 3 0.1 7.1 7.4
```

The program:
1. **computes** d and times the loop;
2. **prints** the first 5 elements;
3. **Check 1:** verifies that *every* `d[i]` equals `a*x + y` as the program computes it;
4. **Check 2** (if `expected` is given): compares every `d[i]` with `expected` in two ways:
   - **exactly**, with `==`;
   - **with a tolerance**: `|d - expected| <= 4·DBL_EPSILON·|expected|`.

It exits with code `0` if all checks pass (within tolerance), `1` on wrong input or not enough
memory, and `2` if a check fails.

## Results (a = 3, x = 0.1, y = 7.1, expected = 7.4)

| N | Memory (3 vectors) | Time of the loop | Check 1 | `d == 7.4` | Within tolerance |
|---|---|---|---|---|---|
| 10 | ~0 MB | ~0 s | ✅ | ❌ | ✅ |
| 10⁶ | 24 MB | ~3 ms | ✅ | ❌ | ✅ |
| 10⁸ | 2.4 GB | ~0.4 s | ✅ | ❌ | ✅ |

### Did you find any problems running the code for some N?
Yes, at **N = 10⁸**. The three vectors of `double` (8 bytes each) need
3 × 10⁸ × 8 B = **2.4 GB of RAM**.
- On a machine or VM with little memory, `malloc` fails, or the system kills the process.
  The program now checks the result of `malloc` and stops with a clear error message.
- Printing every element (as the first version did) writes 10⁸ lines, about 2 GB of text, and
  takes far longer than the computation itself. Version 2 prints only the first elements and
  checks the rest automatically.

### Were you able to test the sum and product of point 4 correctly?
**Not with `==`.** In binary floating point, `3*0.1 + 7.1` gives `7.3999999999999995`,
while the literal `7.4` is stored as `7.4000000000000004`.
- 0.1, 7.1 and 7.4 have no exact binary representation, just like 1/3 has no exact decimal one.
- Every operation rounds the result to the nearest representable `double`.

The two values differ by 8.9 × 10⁻¹⁶, the smallest possible step between doubles near 7.4.
The test therefore has to allow for rounding: compare with a **relative tolerance** of a few
machine epsilons (`DBL_EPSILON` ≈ 2.2 × 10⁻¹⁶). With that tolerance, all elements pass.
