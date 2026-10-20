# HomeWork 03: Matrix multiplication

`MatrixMult.c` computes **C = A B**, where

  C<sub>ij</sub> = Σ<sub>k</sub> A<sub>ik</sub> B<sub>kj</sub>   (i, j, k = 0 … N−1),

with A = a·𝟙<sub>N</sub> and B = b·𝟙<sub>N</sub>.

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
| `N` | Matrix dimension (N × N) |
| `fileout` | File where C is saved (one row per line) |
| `ones` *(default)* | 𝟙<sub>N</sub> = matrix of ones: every A<sub>ij</sub> = a, every B<sub>ij</sub> = b, so **every C<sub>ij</sub> = N·a·b** |
| `identity` | 𝟙<sub>N</sub> = identity matrix: A = a·I, B = b·I, so **C = a·b·I** |

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
`==`: each C<sub>ij</sub> is a sum of N terms, and each addition can add a rounding error, so
the total error can reach ≈ N·ε·|C<sub>ij</sub>| (ε = `DBL_EPSILON` ≈ 2.2·10⁻¹⁶).

For example, `./matrixmult 0.1 0.7 1000 C.dat` gives every C<sub>ij</sub> = 69.999999999999659
instead of 70:
- with `==`, **0 of the 10⁶** elements are exactly equal to 70;
- with the tolerance, all of them pass.

**Is there a quicker way, without comparing every element?** Yes: compare **one number
instead of N²**. The sum of all the elements of C is known in advance:
- `ones`: Σ<sub>ij</sub> C<sub>ij</sub> = N² · N·a·b = **N³·a·b**
- `identity`: Σ<sub>ij</sub> C<sub>ij</sub> = trace = **N·a·b**

So one sum and **a single comparison** are enough.

The limitation: errors that cancel out (one element +1 too big, another −1 too small) go
unnoticed. In any case, checking costs O(N²) and the product costs O(N³), so even the full
check is negligible for large N.

> Another observation: in the `ones` case, every C<sub>ij</sub> is computed with *exactly the
> same* sequence of floating-point operations, so all elements are bit-for-bit identical.
> In principle, checking C<sub>00</sub> plus the fact that all elements are equal is enough.

## Benchmark: which loop ordering is the fastest and why?

The three loops `i, j, k` can be nested in 6 orders. The arithmetic is identical; only the
**order in which memory is visited** changes.

Results with N = 1000 (`-O3`, this PC):

| Ordering | Time [s] | GFLOP/s | Innermost loop walks along… |
|---|---|---|---|
| **ikj** | **0.39** | **5.1** | a row of B and a row of C: contiguous O |
| **kij** | 0.42 | 4.7 | a row of B and a row of C: contiguous O |
| jik | 1.08 | 1.9 | a row of A (contiguous) and a column of B (strided) |
| ijk | 1.19 | 1.7 | a row of A (contiguous) and a column of B (strided) |
| kji | 7.69 | 0.26 | a column of A and a column of C: strided X |
| jki | 7.73 | 0.26 | a column of A and a column of C: strided X |

**ikj / kij are fastest, about 20× faster than jki / kji.**

**Why?**
- **C stores matrices row by row** (row-major). Element (i, j) is at position `i*N + j`, so
  moving along **j** means moving to the next memory address. Moving along **i** or **k** in
  the first index means jumping N·8 bytes = 8 kB.
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
- Each ordering starts from C = 0 and is checked afterwards, so a "fast but wrong" ordering
  would be caught.
- The numbers vary slightly between runs. For a careful benchmark, repeat each measurement
  several times and take the minimum or the mean.
