# HomeWork 03: Matrix multiplication

`MatrixMult.c` computes $C = AB$, where

$$C_{ij} = \sum_k A_{ik} B_{kj}, \qquad i, j, k = 0, \dots, N-1,$$

with $A = a \mathbb{1}_N$ and $B = b \mathbb{1}_N$.

The program:
- computes the product with **all 6 loop orderings** and benchmarks them;
- checks the result;
- saves C to a file.

| File | Content |
|---|---|
| `MatrixMult.c` | Solution of the homework |
| `Makefile` | Builds `matrixmult` and prints how to use it (`make help` prints it again) |

## Build and run

```bash
make
./matrixmult a b N fileout [ones|identity]
```

| Argument | Meaning |
|---|---|
| `a`, `b` | Value of the elements of A and B |
| `N` | Matrix dimension ($N \times N$) |
| `fileout` | File where C is saved (one row per line) |
| `ones` *(default)* | $\mathbb{1}_N$ = matrix of ones: every $A_{ij} = a$, every $B_{ij} = b$, so every $C_{ij} = N a b$ |
| `identity` | $\mathbb{1}_N$ = identity matrix: $A = a I$, $B = b I$, so $C = a b I$ |

Examples:
```bash
./matrixmult 2 3 4 C.dat          # small: also prints A, B and C on screen
./matrixmult 2 3 1000 C.dat       # benchmark (~20 s in total)
make run                          # same as the line above
```

Exit code: `0` all checks passed, `1` wrong input / no memory / file not writable,
`2` a check failed.

## Checking the result

**Full check.** Every element is compared with the expected value. With a tolerance, not
`==`: each $C_{ij}$ is a sum of $N$ terms, and each addition can add a rounding error, so
the total error can reach $\approx N \varepsilon \lvert C_{ij} \rvert$ ($\varepsilon$ = `DBL_EPSILON` ≈ 2.2·10⁻¹⁶).

For example, `./matrixmult 0.1 0.7 1000 C.dat` gives every $C_{ij}$ = 69.999999999999659
instead of 70:
- with `==`, **0 of the 10⁶** elements are exactly equal to 70;
- with the tolerance, all of them pass.

**Is there a quicker way, without comparing every element?** Yes: compare **one number
instead of $N^2$**. The sum of all the elements of C is known in advance:
- `ones`: $\sum_{ij} C_{ij} = N^2 \cdot N a b = N^3 a b$
- `identity`: $\sum_{ij} C_{ij} = \operatorname{tr} C = N a b$

So one sum and **a single comparison** are enough.

The limitation: errors that cancel out (one element $+1$ too big, another $-1$ too small) go
unnoticed. In any case, checking costs $O(N^2)$ and the product costs $O(N^3)$, so even the full
check is negligible for large $N$.

> Another observation: in the `ones` case, every $C_{ij}$ is computed with *exactly the
> same* sequence of floating-point operations, so all elements are bit-for-bit identical.
> In principle, checking $C_{00}$ plus the fact that all elements are equal is enough.

## Benchmark: which loop ordering is the fastest and why?

The three loops `i, j, k` can be nested in 6 orders. The arithmetic is identical; only the
**order in which memory is visited** changes.

**Test conditions**

| | |
|---|---|
| Command | `./matrixmult 2 3 1000 C.dat`, i.e. A and B are 1000 × 1000 matrices with every element equal to a = 2 and b = 3 |
| Data | `double` (8 bytes): 3 matrices × 1000 × 1000 × 8 bytes = 24 MB in total |
| Compiler | gcc 11.4, flags `-Wall -Wextra -O3` (as in the `Makefile`) |
| Machine | laptop with an Intel Core i7-13700H (24 MB of L3 cache), Ubuntu 22.04 under WSL2 (Windows Subsystem for Linux), one core (the code is serial) |
| Timing | `clock_gettime(CLOCK_MONOTONIC)` around the multiplication only; the program was run 5 times and the table shows the **best** of the 5 times for each ordering |
| Work | each product is $2N^3 = 2 \times 10^9$ floating-point operations (one multiplication and one addition per term), so GFLOP/s = 2 / time [s] |

The absolute times depend on the machine and change a little from run to run (between the
5 runs by up to about 20%); the **ranking and the ratios** between orderings are stable.

**Results**

| Ordering (outer → inner loop) | Time [s] | GFLOP/s | What the innermost loop reads | Memory access |
|---|---|---|---|---|
| **ikj** | **0.148** | **13.5** | a row of B and a row of C | $\color{green}{\textbf{contiguous}}$ |
| **kij** | 0.169 | 11.8 | a row of B and a row of C | $\color{green}{\textbf{contiguous}}$ |
| jik | 0.413 | 4.8 | a row of A and a column of B | half contiguous, half strided |
| ijk | 0.486 | 4.1 | a row of A and a column of B | half contiguous, half strided |
| kji | 2.80 | 0.71 | a column of A and a column of C | $\color{red}{\textbf{strided}}$ |
| jki | 2.86 | 0.70 | a column of A and a column of C | $\color{red}{\textbf{strided}}$ |

Meaning of the last column:
- **contiguous**: the innermost loop moves to the *next* memory address at every step, so each
  cache line loaded is fully used (the fast case).
- **strided**: the innermost loop *jumps* $N = 1000$ elements (8 kB) at every step, so each
  access needs a new cache line (the slow case).
- **half contiguous, half strided**: one of the two matrices is read contiguously, the other with
  jumps.

**ikj and kij are the fastest: about 20× faster than jki / kji and about 3× faster than
ijk / jik, although all six do exactly the same arithmetic.**

**Why?**
- **C stores matrices row by row** (row-major). Element (i, j) is at position `i*N + j`, so
  moving along **j** means moving to the next memory address. Moving along **i** or **k** in
  the first index means jumping $8N$ bytes = 8 kB.
- **The CPU reads memory in cache lines** of 64 bytes, i.e. 8 doubles at a time:
  - When the innermost loop moves along j (ikj, kij), each cache line loaded is fully used, and
    the hardware prefetcher can predict the next one.
  - The compiler can also vectorize the loop (SIMD: several multiply-adds per instruction).
- **In jki / kji the innermost loop moves along i.** Every access to A and C jumps to a
  different row, so it is a new cache line from which only 1 of the 8 doubles is used. Most of
  the time is spent waiting for memory, not computing.
- **ijk / jik are in between.** A is read contiguously, but B is read down a column.

**Rule of thumb:** with row-major arrays, put the index that runs along the rows (the last
index, here `j`) in the **innermost** loop.

**How it is timed.**
- `clock_gettime(CLOCK_MONOTONIC)` measures the real elapsed time around the multiplication
  only. It doesn't include filling the matrices, the checks or the file output.
- Each ordering starts from $C = 0$ and is checked afterwards, so a "fast but wrong" ordering
  would be caught.
- The numbers vary between runs (other programs, CPU frequency changes). That's why the program
  was run 5 times and the best time is reported: the minimum is the run least disturbed by
  the rest of the system.
