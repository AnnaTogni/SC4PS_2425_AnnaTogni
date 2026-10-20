# HomeWork 06: Split the work

| Folder | Content |
|---|---|
| [`chunked_daxpy/`](chunked_daxpy/README.md) | d = a·x + y computed in chunks, with partial sums, input from a `Variable = Value` file (GLib), output in HDF5, and 24 tests (`make test`) |

**Main results:**
- The chunked d is **bit-for-bit identical** to the single-loop d.
- The sum of the partial chunk sums is **not** always identical to the single-loop sum, because
  the numbers are added in a different order. It is usually **more accurate**: for n = 10⁷ and
  chunks of 1000, the error drops from 3·10⁻⁴ to 1·10⁻⁹.
