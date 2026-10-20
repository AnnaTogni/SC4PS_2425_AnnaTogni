# Errors found in myfft_tutorial, and proposed solutions

Project: <https://github.com/alexninogh/myfft_tutorial> (commit `ae4095d`).

Each section below is one GitHub issue: the title, followed by the text of the issue.
All problems were reproduced, and all fixes were tested. The complete fix for issues 1–3 is in
[`fixes.patch`](fixes.patch).

| # | Type | Title |
|---|---|---|
| 1 | Bug | Out-of-bounds read in `report_mode_near_frequency()` when a normal mode is above Nyquist |
| 2 | Bug | Python spectrum plot hard-codes the spring constants |
| 3 | Plots | Time-domain sampling plots don't show aliasing (`every 2` drops half the samples) |
| 4 | Docs | Small inaccuracies in README / LESSONS |
| 5 | Code quality | Minor improvements |

---

## Issue 1: Out-of-bounds read in `report_mode_near_frequency()` when a normal mode is above Nyquist

**Description**

In `src/coupled_oscillators_fft.c`, `report_mode_near_frequency()` computes the bin of the
theoretical frequency without checking that it lies inside the spectrum:

```c
best_bin = (size_t)llround(target_frequency_hz / (frequencies[1] - frequencies[0]));
...
if (amplitudes[i] > amplitudes[best_bin])   /* best_bin can be >= bins */
```

With `dt = 0.125 s` the Nyquist frequency is 4 Hz (`bins = 513`). If a mode is above 4 Hz,
`best_bin` is larger than the array, and `amplitudes[best_bin]` / `frequencies[best_bin]`
read outside the allocated memory. The worksheet's Part 7 Option C ("change the spring
constants") leads students directly into this case.

**How to reproduce**

Set `.coupling_k = 400.0` (so f_out = √825/2π = 4.57 Hz > 4 Hz) and run:

```
out-of-phase mode: expected 4.57138 Hz, measured 0.00022 Hz, amplitude 0.0000
```

The "measured" value is garbage read from outside the array. Compiling with
`-fsanitize=address` reports:

```
ERROR: AddressSanitizer: heap-buffer-overflow
READ of size 8 ... in report_mode_near_frequency src/coupled_oscillators_fft.c:156
```

**Proposed solution**

Check the target against the Nyquist frequency before using it as an index. This also tells
the student *why* the mode can't be measured, which is a nice extra lesson on aliasing:

```c
const double nyquist_hz = frequencies[bins - 1U];
...
if (target_frequency_hz >= nyquist_hz) {
    printf("  %s: expected %.5f Hz is above the Nyquist frequency %.5f Hz:"
           " it cannot be measured (aliased), reduce dt\n",
           label, target_frequency_hz, nyquist_hz);
    return;
}
```

With the fix: no AddressSanitizer error, and the output is

```
out-of-phase mode: expected 4.57138 Hz is above the Nyquist frequency 4.00000 Hz: it cannot be measured (aliased), reduce dt
```

---

## Issue 2: Python spectrum plot hard-codes the spring constants

**Description**

`scripts/plot_coupled_oscillators.py` redefines the physical parameters instead of reading them:

```python
mass = 1.0
wall_k = 25.0
coupling_k = 7.0
```

When a student follows LESSONS.md Part 7 Option C and changes `wall_k` or `coupling_k` in
`src/coupled_oscillators_fft.c`, the FFT peaks move, but the "theory" lines in
`coupled_oscillators_spectrum.png` stay at the old values. The plot then suggests that theory
and simulation disagree, which is exactly the opposite of the lesson. Parameters should be
defined in one place only.

**How to reproduce**

Set `.coupling_k = 20.0` and run `make plot-python`. The out-of-phase peak is at 1.28 Hz, but
the green "theory" line is still drawn at 0.9939 Hz.

**Proposed solution**

Let the C program, which knows the parameters actually used, write the theoretical
frequencies to `output/coupled_oscillators_modes.csv`:

```
mode,frequency_hz
in_phase,0.795774715459
out_of_phase,0.993922301044
```

Then read that file in the script, instead of recomputing the frequencies:

```python
with (OUTPUT_DIR / "coupled_oscillators_modes.csv").open(newline="", encoding="utf-8") as handle:
    modes = {row["mode"]: float(row["frequency_hz"]) for row in csv.DictReader(handle)}
in_phase_hz = modes["in_phase"]
out_of_phase_hz = modes["out_of_phase"]
```

Tested: with `coupling_k = 20`, the green line moves to 1.2831 Hz and matches the peak.

---

## Issue 3: Time-domain sampling plots don't show aliasing (`every 2` drops half the samples)

**Description**

1. In `scripts/plot_sampling.gp`, the "good sampling" signal is plotted with `every 2`:

   ```gnuplot
   plot "output/good_sampling_signal.csv" using 1:2 every 2 with lines ...
   ```

   This draws only every second sample, i.e. an effective 256 Hz instead of 512 Hz. The plot
   meant to show *good* sampling is down-sampled. Only about 2 points remain per period of
   the 120 Hz component, which is barely above the Nyquist limit.
2. In both the gnuplot and the Python versions, the good and undersampled signals are plotted
   over the full 1 s: 50–120 oscillations in one panel. Both look like an unreadable zigzag, and
   the false 8 Hz oscillation of the undersampled case can't be seen.
3. `make plot` and `make plot-python` write the same files (`plots/*.png`), so each one
   silently overwrites the other's figures.

**Proposed solution**

- Remove `every 2`.
- Plot the samples as points joined by lines (`linespoints`), so the sampling is visible.
- Zoom the time axis: the first 0.1 s for the good sampling, and 0.5 s for the undersampled
  case (a few periods of the 8 Hz alias).

```gnuplot
set xrange [0:0.1]
plot "output/good_sampling_signal.csv" using 1:2 with linespoints lw 2 pt 7 ps 0.5 title "samples"

set xrange [0:0.5]
plot "output/undersampled_signal.csv" using 1:2 with linespoints lw 2 pt 7 ps 0.5 title "samples"
```

Tested: the good-sampling panel shows all the samples, and the undersampled panel clearly
shows the slow false oscillation.

For point 3, write to `plots/gnuplot/` and `plots/python/` (or add a suffix to the file
names).

---

## Issue 4: Small inaccuracies in README / LESSONS

1. **README, "Suggested modifications" no. 3:** *"Remove the mean-subtraction step in
   `src/fft_utils.c` and inspect how the DC component changes the spectra."*

   All the generated signals already have (almost) zero mean:

   | Signal | Mean |
   |---|---|
   | good_sampling | 1·10⁻¹⁴ |
   | undersampled | 4·10⁻¹⁹ |
   | coupled oscillators | 3·10⁻⁵ |
   | short_record | 0.03 |

   So removing the step shows practically no difference. Suggest also adding a constant
   offset to a signal (e.g. `+ 0.5`) to make the DC bin visible.
2. **README, "Requirements":** *"In this environment `gsl-config` is available, while
   `pkg-config gsl` is not configured"* describes the machine where the project was generated,
   not the reader's machine.
3. **LESSONS.md:** some numbered items are statements, not questions: Part 3 items 1 and 4,
   Part 4 item 1, and Option C item 1. SOLUTIONS.md answers them with "Given." They could be
   written as introductory text instead of numbered questions.
4. **Lesson strings:** the lesson strings in `sampling_demo.c` are fixed text. After Part 7
   Option B (N = 512), the program finds both 50 Hz and 55 Hz, but still prints "the record is
   too short, so the FFT cannot separate nearby frequencies".

---

## Issue 5: Minor code-quality improvements

- **Use GSL's real-input FFT.** `fft_utils.c` copies real data into a complex array and calls
  `gsl_fft_complex_radix2_forward`. The real-input transforms (`gsl_fft_real_radix2_transform`,
  or the mixed-radix `gsl_fft_real_transform`) need half the memory and work.
  `gsl_fft_real_transform` also works for any N, not only powers of two. If the power-of-two
  restriction is kept on purpose for teaching, the README should say it's a choice, not a
  requirement of FFTs.
- **Redundant check.** In `sampling_demo.c`, `main()` checks
  `fftu_is_power_of_two(cases[0].sample_count)` before the loop that checks every case anyway.
- **Duplicated code.** `ensure_directory()` and `write_spectrum_csv()` are copied in both
  programs. They could move to a shared utility file, like `fft_utils.c`.
- **Unchecked writes.** `write_*_csv()` don't check the return value of `fclose()` (write
  errors, e.g. a full disk, go unnoticed). `snprintf()` path truncation isn't checked either.
