# Instructor Key

This file contains concise model answers for the guided worksheet in
`LESSONS.md`. The answers are intentionally short. They are meant to support
discussion, not replace it.

## Part 1: Before running anything

1. An FFT estimates how strongly different frequencies are present in a sampled
   signal.
2. Sampling rate is samples per second. Nyquist frequency is half the sampling
   rate and is the highest frequency that can be represented without aliasing:
   $f_N = \dfrac{f_s}{2}$, where $f_s$ is the sampling frequency and $f_N$ is
   the Nyquist limit. Frequency resolution is the spacing between FFT bins:
   $\Delta f = \dfrac{1}{T} = \dfrac{f_s}{N}$, where $\Delta f$ is the bin
   spacing, $T$ is the total acquisition time, and $N$ is the number of
   samples.
3. The FFT routine used here assumes equally spaced samples. If the time grid is
   not uniform, the bin frequencies no longer match the data correctly.

## Part 2: Build and run

1. The generated CSV files are:
   - `output/good_sampling_signal.csv`
   - `output/good_sampling_spectrum.csv`
   - `output/undersampled_signal.csv`
   - `output/undersampled_spectrum.csv`
   - `output/short_record_signal.csv`
   - `output/short_record_spectrum.csv`
   - `output/coupled_oscillators_time.csv`
   - `output/coupled_oscillators_spectrum.csv`
2. `src/sampling_demo.c`
3. `src/coupled_oscillators_fft.c`
4. About `50 Hz` and `120 Hz`

## Part 3: Sampling and aliasing

1. Given.
2. $f_N = \dfrac{f_s}{2} = \dfrac{128\,\mathrm{Hz}}{2} = 64\,\mathrm{Hz}$
3. Because $120\,\mathrm{Hz} > 64\,\mathrm{Hz}$, so it lies above Nyquist and folds into a false
   lower frequency after sampling.
4. Given.
5. A sampled sinusoid is only observed at discrete times. Different continuous
   frequencies can produce the same sampled sequence if they differ by integer
   multiples of the sampling frequency. Here the alias is
   $f_{\mathrm{alias}} = |120 - 128|\,\mathrm{Hz} = 8\,\mathrm{Hz}$, where
   $f_{\mathrm{alias}}$ is the observed aliased frequency.

## Part 4: Frequency resolution

1. Given.
2. $T = \dfrac{N}{f_s} = \dfrac{64}{512}\,\mathrm{s} = 0.125\,\mathrm{s}$,
   where $T$ is the total record length.
3. $\Delta f = \dfrac{1}{T} = 8\,\mathrm{Hz}$, where $\Delta f$ is the FFT bin spacing.
4. The two frequencies differ by only $5\,\mathrm{Hz}$, which is smaller than the FFT bin
   spacing, so they are not cleanly resolved.
5. Change the total acquisition time. Plotting tools and file formats do not
   improve the underlying frequency resolution.

## Part 5: Coupled oscillators

1. $m \ddot{x}_1 = -(k + k_c)x_1 + k_c x_2$ and
   $m \ddot{x}_2 = k_c x_1 - (k + k_c)x_2$, where $m$ is the mass, $k$ is the
   wall-spring constant, $k_c$ is the coupling constant, and $x_1$, $x_2$ are
   the displacements.
2. Two identical masses attached to walls and coupled by a spring
3. An in-phase mode and an out-of-phase mode
4. The chosen initial condition excites a combination of normal modes rather
   than only one of them.
5. It checks whether the numerical solution is consistent with the analytical
   physics and helps detect modeling or coding mistakes.

## Part 6: Plot inspection

1. The sampling spectra plot, especially the undersampled spectrum
2. The short-record spectrum
3. Yes, the peaks line up closely with the theoretical markers.
4. Small differences come from finite observation time, discrete FFT bin
   spacing, solver tolerance, and the fact that a peak rarely lands exactly on a
   single bin.

## Part 7: Small code modifications

### Option A: Change the sampling rate

Model answer:

1. Aliasing occurs if a true frequency exceeds `fs / 2`.
2. Predicted peak positions should be compared against the Nyquist limit
   $f_N = \dfrac{f_s}{2}$ and possible fold back by multiples of $f_s$, where
   $f_s$ is the sampling frequency.
3. The FFT output should match the prediction within the FFT bin spacing.

### Option B: Change the observation time

Model answer:

1. Increasing the number of samples while keeping $\Delta t$ fixed increases
   $T = N \Delta t$.
2. Since $\Delta f = \dfrac{1}{T}$, the bin spacing becomes smaller.
3. The $50\,\mathrm{Hz}$ and $55\,\mathrm{Hz}$ peaks should become easier to distinguish.

### Option C: Change the spring constants

Model answer:

1. Given.
2. The in-phase mode scales like
   $f_{\mathrm{in}} = \dfrac{1}{2\pi}\sqrt{\dfrac{k}{m}}$, while the
   out-of-phase mode scales like
   $f_{\mathrm{out}} = \dfrac{1}{2\pi}\sqrt{\dfrac{k + 2k_c}{m}}$.
3. The new FFT peaks should shift consistently with those formulas.

## Part 8: Reflection

Accept any reasonable response that includes:

1. one correct FFT concept
2. one correct point about reproducible project structure
3. one plausible improvement or extension

## Optional extension

Any good answer should include:

1. a clearly stated model
2. chosen parameters
3. a physically sensible spectral prediction
4. a comparison between prediction and computed FFT

Examples:

- A damped oscillator should show a broadened peak compared with the undamped
  case.
- A driven oscillator should show a strong component near the drive frequency.
- Added random noise should raise the spectral background and make peaks less
  clean.
