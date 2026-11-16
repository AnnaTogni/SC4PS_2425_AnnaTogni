# Errors found in myfft_tutorial, and proposed solutions

Project: <https://github.com/alexninogh/myfft_tutorial> (commit `ae4095d`).

Each section is one GitHub issue. All problems were reproduced and all fixes were tested.

| # | Type | Title |
|---|---|---|
| 1 | Bug | Out-of-bounds read when a normal mode is above the Nyquist frequency |
| 2 | Bug | The Python spectrum plot hard-codes the spring constants |
| 3 | Plots | The time-domain sampling plots don't show aliasing |
| 4 | Code quality | Minor improvements |

---

## Issue 1: Out-of-bounds read when a normal mode is above the Nyquist frequency

**Problem.** In `src/coupled_oscillators_fft.c`, `report_mode_near_frequency()` turns the
theoretical frequency into a bin index without checking that it is inside the spectrum:

```c
best_bin = (size_t)llround(target_frequency_hz / (frequencies[1] - frequencies[0]));
```

With `dt = 0.125 s` the Nyquist frequency is 4 Hz. If a mode is above 4 Hz, `best_bin` is
past the end of the array and the program reads outside the allocated memory. The worksheet's
Part 7 Option C ("change the spring constants") leads straight into this case.

**How to reproduce.** Set `.coupling_k = 400.0` (out-of-phase mode at 4.57 Hz) and run: the
"measured" frequency is a random value, and compiling with `-fsanitize=address` reports a
`heap-buffer-overflow` in `report_mode_near_frequency`.

**Proposed solution.** Compare the frequency with the Nyquist frequency first, and explain why
it can't be measured:

```c
const double nyquist_hz = frequencies[bins - 1U];
if (target_frequency_hz >= nyquist_hz) {
    printf("  %s: expected %.5f Hz is above the Nyquist frequency %.5f Hz:"
           " it cannot be measured (aliased), reduce dt\n",
           label, target_frequency_hz, nyquist_hz);
    return;
}
```

Tested: no more memory error, and the program prints the explanation.

---

## Issue 2: The Python spectrum plot hard-codes the spring constants

**Problem.** `scripts/plot_coupled_oscillators.py` repeats the physical parameters
(`mass = 1.0`, `wall_k = 25.0`, `coupling_k = 7.0`) instead of reading them from the C
program. If a student changes the spring constants in the C code (Part 7 Option C), the FFT
peaks move but the "theory" lines in the plot stay where they were. Theory and simulation then
seem to disagree, which is the opposite of the lesson.

**How to reproduce.** Set `.coupling_k = 20.0` and run `make plot-python`: the peak is at
1.28 Hz, but the theory line is still at 0.99 Hz.

**Proposed solution.** Define the parameters in one place only. The C program writes the
theoretical frequencies to `output/coupled_oscillators_modes.csv`, and the script reads them:

```python
with (OUTPUT_DIR / "coupled_oscillators_modes.csv").open(newline="", encoding="utf-8") as handle:
    modes = {row["mode"]: float(row["frequency_hz"]) for row in csv.DictReader(handle)}
```

Tested: with `coupling_k = 20` the theory line moves to 1.28 Hz and matches the peak.

---

## Issue 3: The time-domain sampling plots don't show aliasing

**Problem.**
1. `scripts/plot_sampling.gp` plots the "good sampling" signal with `every 2`, i.e. only every
   second sample: the plot meant to show *good* sampling is itself undersampled.
2. Both signals are drawn over the full 1 s (50 to 120 oscillations in one panel), so they look
   like an unreadable zigzag and the false 8 Hz oscillation can't be seen.
3. `make plot` and `make plot-python` write to the same files, so each one overwrites the
   other's figures.

**Proposed solution.** Remove `every 2`, draw the samples as points (`linespoints`) and zoom in
on the time axis (0.1 s for the good sampling, 0.5 s for the undersampled case):

```gnuplot
set xrange [0:0.1]
plot "output/good_sampling_signal.csv" using 1:2 with linespoints pt 7 ps 0.5 title "samples"
```

Write the two sets of plots to different folders (`plots/gnuplot/`, `plots/python/`).

Tested: the undersampled panel now clearly shows the slow false oscillation.

---

## Issue 4: Minor code-quality improvements

- **Real-input FFT.** `fft_utils.c` copies real data into a complex array. GSL's real-input
  FFT (`gsl_fft_real_transform`) needs half the memory and also works when N is not a power
  of two.
- **Duplicated code.** `ensure_directory()` and `write_spectrum_csv()` are copied in both
  programs; they could go in a shared file such as `fft_utils.c`.
- **Unchecked writes.** The return value of `fclose()` is not checked, so a write error (e.g. a
  full disk) goes unnoticed.
