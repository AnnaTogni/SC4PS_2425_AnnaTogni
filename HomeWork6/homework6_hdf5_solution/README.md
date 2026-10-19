# HomeWork 06 HDF5 Solution

This directory contains the C + HDF5 solution for the assignment in `HomeWork06.pdf`.

## What the program does

- Computes the original `d[i] = a * x[i] + y[i]` with one loop.
- Recomputes the same result by splitting the vector into chunks.
- Builds `partial_chunk_sum` by summing each chunk separately.
- Verifies that the chunked result matches the original vector and total sum.
- Prints the execution time of the reference loop, the chunked loop, and the HDF5 write step.
- Supports the same interactive input style as `NumericalSum.c`.
- Still supports the required plain-text `Variable = Value` configuration file mode.
- Writes vectors, chunk datasets, and partial sums to an HDF5 file when HDF5 is available.

## Files

- `main.c`: the homework implementation.
- `homework6.cfg`: sample input file.
- `Makefile`: build and run helper.
- `results.h5`: example HDF5 output file.
- `explanation.html`: step-by-step explanation of the code for a physics PhD student.

## Build

```bash
cd homework6_hdf5_solution
make
```

## Run

```bash
./homework6
```

This matches `NumericalSum.c`:

- If you launch `./homework6 N a`, the program uses those two values directly.
- If you launch `./homework6`, it asks `Enter vector size N:` and `Enter constant a:`.
- Then it asks the same `x` and `y` questions as `NumericalSum.c`, plus the extra homework inputs `chunk size` and `HDF5 output file name`.
- For the HDF5 file name, you can type `out` or `out.h5`. If you omit the extension, the program now adds `.h5` automatically.

## Config-file mode

```bash
./homework6 --config homework6.cfg
```

## To read .h5 files

h5ls -r : shows the tree of groups and datasets
h5dump -n : shows names only
h5dump : prints the actual stored values

```bash
h5ls -r output_chunks.h5
h5dump -n output_chunks.h5
h5dump output_chunks.h5
```


## HTML explanation

Open `explanation.html` in a browser to read a guided explanation of the code structure,
the numerical logic, the chunking strategy, and the HDF5 output.

## Notes

- Configuration parsing uses `glib`'s `GKeyFile`. The code wraps the assignment-style `key = value` file with a temporary section header internally so the input file can stay in the required format.
- In config mode, vectors are provided explicitly as `x_0`, `x_1`, ..., `y_0`, `y_1`, ...
- The program calculates `number_of_chunks` automatically from `vector_length` and `chunk_size`.
- HDF5 output is enabled automatically when HDF5 headers and libraries are installed. On Ubuntu this usually means installing `libhdf5-dev`.
