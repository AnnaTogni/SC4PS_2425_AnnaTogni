# HomeWork 06: Split the work

The daxpy code of HomeWork 2, **d = a·x + y**, with the same vectors: a is a scalar, all the
elements of x are equal to x and all the elements of y are equal to y (here a = 3, x = 0.1,
y = 7.1, as in HomeWork 2). It is rewritten in C to follow the four points of the assignment:

1. compute d **chunk by chunk** instead of with a single loop, and check that d is the same;
2. sum the elements of every chunk into `partial_chunk_sum[]`, sum the partial sums, and check
   that the result is the same as the sum of d in the original code;
3. take **all the inputs from an input file**, a text file with one `Variable = Value` line per
   input, read with an existing configuration-file library;
4. save **all the vector chunks and the partial sums in an HDF5 file**.

| File | Content |
|---|---|
| `chunked_daxpy.c` | The program |
| `daxpy.cfg` | The input file |
| `Makefile` | `make` (builds and prints the instructions), `make run`, `make help`, `make clean` |

Requirements: `gcc`, GLib (`libglib2.0-dev`, for the configuration-file parser), HDF5
(`libhdf5-dev`).

## Build and run

```bash
make
./chunked_daxpy daxpy.cfg      # or simply ./chunked_daxpy (daxpy.cfg is the default)
h5ls -r daxpy_chunks.h5        # what is inside the HDF5 file
h5dump -d /chunks/d daxpy_chunks.h5
```

Output with the given `daxpy.cfg`:

```
n = 100, chunk_size = 8, number_of_chunks = 13 (last chunk: 4 elements)
a = 3, every x[i] = 0.1, every y[i] = 7.1

1. d (chunks) == d (single loop), element by element: YES (0 of 100 differ)

2. sum of d
   single loop:                 739.99999999999864
   sum of partial_chunk_sum[]:  740.00000000000011
   exactly equal: NO;  equal within rounding (|diff| <= n eps sum|d| = 1.6e-11): YES

4. HDF5 file 'daxpy_chunks.h5' written (inspect: h5ls -r daxpy_chunks.h5)
```

## 1. Split the work in chunks

```c
number_of_chunks = (n + chunk_size - 1) / chunk_size;      // integer CEILING(n / chunk_size)
for (chunk = 0; chunk < number_of_chunks; chunk++) {
    start = chunk * chunk_size;                            // 0-based
    end   = min(start + chunk_size, n);                    // one past the last element
    for (i = start; i < end; i++)
        d[i] = a * x[i] + y[i];
}
```

The assignment warns that its pseudocode "may contain errors". Two points must be changed when
writing it in C:
1. **CEILING(n / chunk_size):** in C, `n / chunk_size` between two integers is *already
   truncated*, so `ceil(n / chunk_size)` gives 12 instead of 13 for n = 100, chunk_size = 8. The
   integer formula `(n + chunk_size − 1) / chunk_size` gives the correct 13.
2. **current_start = chunk_index · chunk_size + 1:** the pseudocode counts from 1, C arrays from
   0. Copied literally, it skips element 0 and writes one element past the end of the array. In C
   the chunk k covers the elements from `k·chunk_size` to `min((k+1)·chunk_size, n) − 1`; with
   n = 100 and chunk_size = 8 there are 12 chunks of 8 elements and a last one of 4.

**Check: d is the same as in the original code.** The program also computes d with the original
single loop and compares the two vectors element by element. Every element is computed with
exactly the same operation, `a*x[i] + y[i]`, so the two vectors are **bit-for-bit identical**
(0 of 100 elements differ).

## 2. Partial sums

Inside the loop over the chunks, the elements of each chunk are added into
`partial_chunk_sum[chunk]`; then the partial sums are added together.

**Check: is this sum the same as the sum of d in the original code?** Not always bit-for-bit,
and that is correct. The single loop computes ((d₀ + d₁) + d₂) + …, while the chunked version
first adds inside each chunk and then adds the partial sums: these are the **same numbers added
in a different order**, and floating-point addition is not associative (see HomeWork 5). With
the given input file every dᵢ is 3·0.1 + 7.1 = 7.3999999999999995 (as in HomeWork 2), and the two
sums of the 100 elements are 739.99999999999864 and 740.00000000000011: equal except in the last
digits.

So the program checks that they agree **within the rounding error** of a sum of n terms,
|difference| ≤ n·ε·Σ|dᵢ| (here 1.6·10⁻¹¹), instead of comparing them with `==`.

## 3. Input file

All the inputs of the program are read from a text input file, here `daxpy.cfg`, with one
`Variable = Value` line per input. Nothing is written in the code: to run a different case, edit
the file and run the program again, without recompiling. Lines starting with `#` are comments.

```ini
# vector length
n = 100

# give EITHER chunk_size OR number_of_chunks
chunk_size = 8

# d[i] = a * x[i] + y[i]: a is a scalar, every element of x is x, every element of y is y
# (the values of HomeWork 2)
a = 3
x = 0.1
y = 7.1

# HDF5 output file (.h5 is added if missing)
output = daxpy_chunks.h5
```

| Variable | Meaning |
|---|---|
| `n` | length of the vectors |
| `chunk_size` *or* `number_of_chunks` | how to split the work: give one of the two, the other is computed |
| `a` | the scalar in d = a·x + y |
| `x` | the value of every element of the vector x |
| `y` | the value of every element of the vector y |
| `output` | name of the HDF5 file (optional, default `daxpy_chunks.h5`) |

**Reading the file.** As the assignment suggests, the file is read with an existing library,
GLib's `GKeyFile`, which parses `key = value` files. `GKeyFile` expects the keys to be inside a
`[group]` section, so the program adds a `[daxpy]` line in front of the file's content before
parsing it.

**Mistakes in the file are reported, never silently ignored:** a missing variable, a value that
is not a number, a non-positive or non-integer `n` or `chunk_size`, both `chunk_size` and
`number_of_chunks` given (or neither) stop the program with an error message; an unknown
variable name, i.e. a typo, gives a warning.

## 4. HDF5 file

```
/                    attributes: n, chunk_size, number_of_chunks, a
/x  /y  /d           the full vectors (n values)
/partial_chunk_sum   one value per chunk
/total_sum           attributes: chunked (sum of the partial sums), reference (single loop)
/chunks/x            2-D [number_of_chunks][chunk_size]: row k = chunk k of x
/chunks/y            same for y
/chunks/d            same for d   (last row padded with NaN if shorter)
/chunks/start        first element of each chunk
/chunks/length       number of elements of each chunk
```

- The input parameters are stored as **attributes**, so the file says how it was produced.
- Every vector is saved **chunk by chunk** as a 2-D dataset with one row per chunk. When the
  last chunk is shorter, the rest of its row is filled with NaN, and `/chunks/length` gives the
  true length of every chunk.

```
$ h5dump -d /chunks/d -d /chunks/length file.h5      (n = 25, chunk_size = 8)
   (0,0): 3, 3, 3, 3, 3, 3, 3, 3,
   (1,0): 3, 3, 3, 3, 3, 3, 3, 3,
   (2,0): 3, 3, 3, 3, 3, 3, 3, 3,
   (3,0): 3, nan, nan, nan, nan, nan, nan, nan
   length: 8, 8, 8, 1
```

**Why one row per chunk, and not one group per chunk?** Storing each chunk as its own group,
`/chunks/chunk_0000/…`, is fine for 13 chunks, but with millions of chunks it creates millions of
HDF5 objects: the file becomes huge and slow to write. A 2-D dataset stores the same information
as a single object.
