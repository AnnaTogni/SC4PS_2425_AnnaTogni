# Split the work: daxpy in chunks

HomeWork 06: compute **d = a·x + y** chunk by chunk instead of with a single loop:
- sum every chunk into `partial_chunk_sum[]`;
- read every input from a `Variable = Value` file;
- save the chunks and partial sums in an **HDF5** file.

| File | Content |
|---|---|
| `chunked_daxpy.c` | The program |
| `daxpy.cfg` | Example input file |
| `Makefile` | `make` (builds and prints the instructions), `make run`, `make test`, `make help`, `make clean` |
| `tests/run_tests.sh` | 24 tests: edge cases, HDF5 content, wrong input files (`make test`) |

Requirements: `gcc`, GLib (`libglib2.0-dev`, for the config parser), HDF5 (`libhdf5-dev`).

## Build and run

```bash
make
./chunked_daxpy daxpy.cfg      # or simply ./chunked_daxpy (daxpy.cfg is the default)
make test
h5ls -r daxpy_chunks.h5        # what is inside the HDF5 file
h5dump -d /chunks/d daxpy_chunks.h5
```

## 3. Input file

All inputs come from a text file with one `Variable = Value` per line. Lines starting with `#`
are comments.

```ini
n = 100            # vector length (1e6 is accepted too)
chunk_size = 8     # EITHER chunk_size ...
# number_of_chunks = 13   # ... OR number_of_chunks
a = 3
x = 0.1            # x[i] = x + i * x_step
x_step = 0.01      # optional, default 0
y = 7.1            # y[i] = y + i * y_step
y_step = 0         # optional, default 0
output = daxpy_chunks.h5   # optional (.h5 is added if missing)
```

**Parsing.** The file is parsed with GLib's `GKeyFile`, an existing, well-tested library, as
the homework suggests. `GKeyFile` expects `[group]` headers, so the program adds a `[daxpy]`
line in front of the file's content.

**Mistakes in the file are reported, never silently ignored:**
- a missing variable;
- a value that isn't a number;
- a non-positive or non-integer `n` or `chunk_size`;
- both `chunk_size` and `number_of_chunks` given, or neither;
- an unknown variable name, i.e. a typo (this one is a warning).

## 1. Chunks

```c
number_of_chunks = (n + chunk_size - 1) / chunk_size;      // integer CEILING(n / chunk_size)
for (chunk = 0; chunk < number_of_chunks; chunk++) {
    start = chunk * chunk_size;                            // 0-based
    end   = min(start + chunk_size, n);                    // one past the last element
    for (i = start; i < end; i++)
        d[i] = a * x[i] + y[i];
}
```

**The two traps in the homework's pseudocode** (it says "it may contain errors"):
1. **CEILING(n / chunk_size):** in C, `n / chunk_size` between integers is *already truncated*,
   so `ceil(n / chunk_size)` gives 12 instead of 13 for n = 100, chunk_size = 8. The integer
   formula `(n + chunk_size − 1) / chunk_size` is correct.
2. **current_start = chunk_index · chunk_size + 1:** the pseudocode counts from 1, C arrays from
   0. Copied literally, it skips element 0 and writes one element past the end of the array.

**Check 1: d is identical.** Every element of d is computed with exactly the same operations,
`a*x[i] + y[i]`, whether in one loop or in chunks. So the two vectors are **bit-for-bit equal**,
and the program compares them with `==`, not with a tolerance.

## 2. Partial sums

**Check 2: the sums are *not* always identical, and that is correct.** The single loop computes
((d₀ + d₁) + d₂) + …, while the chunked version first adds inside each chunk, then adds the
partial sums. These are the **same numbers added in a different order**, and floating-point
addition is not associative (see HomeWork5). So the two sums can differ in the last digits.

The program therefore checks that they agree **within the rounding error** of a sum of n terms,
|difference| ≤ n·ε·Σ|dᵢ|. It also computes an accurate reference sum (Neumaier compensated
summation) to see which of the two is closer.

Results for n = 10⁷, a = 1, x = 0.1, y = 0.2, i.e. dᵢ = 0.30000000000000004, so the exact
sum is 3000000.0000000005:

| chunk_size | Single-loop sum | Error | Chunked sum | Error |
|---|---|---|---|---|
| 1 | 2999999.9996692175 | −3.3·10⁻⁴ | 2999999.9996692175 | −3.3·10⁻⁴ |
| 8 | 2999999.9996692175 | −3.3·10⁻⁴ | 2999999.9999314556 | −6.9·10⁻⁵ |
| 1000 | 2999999.9996692175 | −3.3·10⁻⁴ | **3000000.0000000019** | **1.4·10⁻⁹** |
| 10⁷ (one chunk) | 2999999.9996692175 | −3.3·10⁻⁴ | 2999999.9996692175 | −3.3·10⁻⁴ |

**The chunked sum is *more accurate*:** with chunk_size = 1000, it is 200 000× more accurate.

- **Single loop:** the running sum quickly becomes large, so each new small term dᵢ is added to a
  big number and loses many of its last bits.
- **Chunks:** each partial sum adds only `chunk_size` terms (small numbers added to small
  numbers), and then only `number_of_chunks` partial sums are added together. This is the idea
  behind *pairwise summation*.
- **The two extremes** (chunk_size = 1, or one single chunk) are just the single loop again.

**Timing** (n = 10⁷, `-O2`): the single loop takes ~0.03 s. The chunked version is about as
fast for large chunks, but slower for chunk_size = 1 (0.13 s), where the loop overhead per
chunk dominates.

## 4. HDF5 file

```
/                    attributes: n, chunk_size, number_of_chunks, a
/x  /y  /d           the full vectors (n values)
/partial_chunk_sum   one value per chunk
/total_sum           attributes: chunked, reference, accurate
/chunks/x            2-D [number_of_chunks][chunk_size]: row k = chunk k of x
/chunks/y            same for y
/chunks/d            same for d   (last row padded with NaN if shorter)
/chunks/start        first element of each chunk
/chunks/length       number of elements of each chunk
```

```
$ h5dump -d /chunks/d -d /chunks/length file.h5      (n = 25, chunk_size = 8)
   (0,0): 3, 3, 3, 3, 3, 3, 3, 3,
   (1,0): 3, 3, 3, 3, 3, 3, 3, 3,
   (2,0): 3, 3, 3, 3, 3, 3, 3, 3,
   (3,0): 3, nan, nan, nan, nan, nan, nan, nan
   length: 8, 8, 8, 1
```

**Why one row per chunk, and not one group per chunk?** A first version stored each chunk as
its own group, `/chunks/chunk_0000/…`. That is fine for 13 chunks, but with 10⁷ chunks it
created a 3.3 GB file and still hadn't finished after several minutes. A 2-D dataset stores
the same information with a single object. With 10⁷ chunks of size 1, the file is 0.7 GB and
is written in 1.3 s.

## Tests (`make test`)

`tests/run_tests.sh` checks, with 24 tests:
- **Correct inputs:**
  - n divisible by chunk_size;
  - not divisible (last chunk shorter);
  - chunk_size > n (one chunk);
  - chunk_size = 1;
  - `number_of_chunks` instead of `chunk_size`;
  - n = 10⁶ with non-constant x;
  - output name without `.h5`.
- **HDF5 content:** the number of chunks, the chunk lengths (8 8 8 1), the NaN padding of the last
  row, and that the partial and total sums are stored.
- **Wrong inputs:**
  - a missing variable;
  - both chunk keys, or neither;
  - chunk_size = 0;
  - n = 2.5;
  - a = abc;
  - more chunks than elements;
  - a file that doesn't exist.

  Each one must exit with code 1 and an explanatory message.
- **Typo warning:** an unknown variable produces a warning.

The tests start from an empty `tests/tmp/` folder every time, so old output files can't make a
test pass. As a check, running them against a program that does nothing makes all 24 fail.
