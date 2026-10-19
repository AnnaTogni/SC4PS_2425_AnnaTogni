# FFT Tutorial Project in C with GNU GSL

This repository is a compact tutorial project for beginner scientific-computing
PhD physics students. It is intentionally small enough to read in one sitting,
but complete enough to demonstrate a realistic workflow:

- write multi-file C code
- compile it with a `Makefile`
- link an external scientific library, GNU GSL
- generate numerical data in CSV format
- inspect that data with FFT-based reasoning
- plot the results with either `gnuplot` or Python

The project uses a simple but important theme from physics and numerical
analysis: a Fourier transform is only as meaningful as the sampling strategy
used to generate the data.

## Learning goals

After working through the files, you should be able to explain:

- why uniformly spaced samples are required by the FFT routines used here
- how the Nyquist frequency limits what can be reconstructed from sampled data
- why a longer acquisition time improves frequency resolution
- how to connect time-domain simulation and frequency-domain analysis
- how to structure a small reproducible scientific-computing project

## Recommended reading order

If you are using this repository as a tutorial, read the files in this order:

1. `README.md`
2. `Makefile`
3. `include/fft_utils.h`
4. `src/fft_utils.c`
5. `src/sampling_demo.c`
6. `src/coupled_oscillators_fft.c`
7. plotting scripts in `scripts/`

This order mirrors a typical scientific workflow: understand the project,
understand the build, understand the reusable numerical utilities, then inspect
the physical applications and visualization layer.

## Repository map

- `Makefile`
  Builds the programs, runs them, and generates plots.
- `LESSONS.md`
  Short worksheet for a lab session or self-study assignment.
- `SOLUTIONS.md`
  Concise instructor key with model answers for the worksheet.
- `include/fft_utils.h`
  Public interface for the shared FFT helper routines.
- `src/fft_utils.c`
  Implementation of the FFT helper routines using GNU GSL.
- `src/sampling_demo.c`
  Teaches sampling, aliasing, and frequency resolution.
- `src/coupled_oscillators_fft.c`
  Solves a two-mass spring system and analyzes its normal modes with an FFT.
- `scripts/plot_sampling.gp`
  `gnuplot` visualization for the sampling tutorial.
- `scripts/plot_coupled_oscillators.gp`
  `gnuplot` visualization for the coupled-oscillator tutorial.
- `scripts/plot_sampling.py`
  Python `matplotlib` version of the sampling plots.
- `scripts/plot_coupled_oscillators.py`
  Python `matplotlib` version of the coupled-oscillator plots.
- `output/`
  Generated CSV data files.
- `plots/`
  Generated PNG figures.

## Requirements

You need:

- a C compiler such as `gcc`
- GNU GSL
- optionally `gnuplot`
- optionally Python 3 with `matplotlib`

This project detects GSL with `gsl-config`:

```bash
gsl-config --version
make
make run
```

In this environment `gsl-config` is available, while `pkg-config gsl` is not
configured, so the `Makefile` uses `gsl-config` directly.

## Quick start

Build the C programs:

```bash
make
```

Generate the CSV output:

```bash
make run
```

Generate figures with `gnuplot`:

```bash
make plot
```

Generate figures with Python:

```bash
make plot-python
```

Clean generated files:

```bash
make clean
```

## The central numerical lesson

An FFT does not magically reveal the truth about a signal. It reveals the
frequency content of the discrete signal you actually sampled. That distinction
matters in research work.

Three questions always matter:

1. Did you sample fast enough?
2. Did you sample long enough?
3. Are your samples uniformly spaced?

This project addresses all three.

### Sampling rate and aliasing

If the sampling frequency is $f_s$ (code variable: `fs` or
`sample_rate_hz`), then the Nyquist frequency is

$$
f_N = \frac{f_s}{2}
$$

where:

- $f_s$ is the sampling frequency
- $f_N$ is the Nyquist frequency

Any component above $f_N$ cannot be represented correctly by the sampled data.
Instead it is folded back into the spectrum and appears as a false low-frequency
signal. This is aliasing.

One useful aliasing rule is

$$
f_{\mathrm{alias}} = \left| f - n f_s \right|
$$

where:

- $f$ is the true continuous frequency
- $f_s$ is the sampling frequency
- $n$ is an integer chosen so that $f_{\mathrm{alias}}$ falls into the visible
  frequency range below Nyquist

In `src/sampling_demo.c`:

- `good_sampling` uses `fs = 512 Hz`, so 50 Hz and 120 Hz are both resolved
- `undersampled` uses `fs = 128 Hz`, so 120 Hz is above Nyquist and aliases to
  about 8 Hz

### Record length and frequency resolution

For a total record length $T = N \Delta t$, the frequency spacing of the FFT
bins is

$$
\Delta f = \frac{f_s}{N} = \frac{1}{T}
$$

where:

- $\Delta f$ is the spacing between FFT frequency bins
- $f_s$ is the sampling frequency
- $N$ is the number of samples
- $\Delta t$ is the sampling time step (code variable: `dt`)
- $T$ is the total observation time, with $T = N \Delta t$

This means that frequency resolution is controlled by the total observation
time, not just the step size. A smaller $\Delta t$ helps extend the frequency
range, but a longer $T$ is what helps distinguish nearby peaks.

In `src/sampling_demo.c`:

- `short_record` uses only 64 samples at `fs = 512 Hz`
- the total duration is only `0.125 s`
- therefore $\Delta f = 8\,\mathrm{Hz}$
- that is too coarse to separate 50 Hz from 55 Hz well

### Uniform output grids matter

The ODE solver in `src/coupled_oscillators_fft.c` uses adaptive internal time
steps for accuracy, but the code still requests the solution on a uniform output
grid before applying the FFT.

This is a realistic scientific-computing pattern:

- use adaptive algorithms internally when they help accuracy
- interpolate or sample the final output onto a uniform grid when the next tool
  requires it

The uniform grid used for the FFT has the form

$$
t_i = i \, \Delta t
$$

where:

- $t_i$ is the $i$-th sample time
- $i$ is the sample index
- $\Delta t$ is the constant sampling step

## Physics example: coupled oscillators

The second tutorial program studies two identical masses connected to walls and
to each other by springs:

$$
m \, \ddot{x}_1 = - (k + k_c) x_1 + k_c x_2
$$

$$
m \, \ddot{x}_2 = k_c x_1 - (k + k_c) x_2
$$

where:

- $m$ is the mass of each oscillator
- $k$ is the spring constant of the wall springs
- $k_c$ is the coupling spring constant
- $x_1(t)$ and $x_2(t)$ are the displacements of the two masses
- $\ddot{x}_1$ and $\ddot{x}_2$ are the second time derivatives, that is,
  accelerations

This system has two normal modes:

- an in-phase mode
- an out-of-phase mode

Their theoretical frequencies are

$$
f_{\mathrm{in}} = \frac{1}{2 \pi} \sqrt{\frac{k}{m}}
$$

$$
f_{\mathrm{out}} = \frac{1}{2 \pi} \sqrt{\frac{k + 2 k_c}{m}}
$$

where:

- $f_{\mathrm{in}}$ is the in-phase normal-mode frequency
- $f_{\mathrm{out}}$ is the out-of-phase normal-mode frequency
- $k$, $k_c$, and $m$ are the same parameters defined above

The code:

1. defines the equations of motion
2. integrates them with a GSL Runge-Kutta solver
3. samples the displacement of one mass uniformly in time
4. computes the FFT amplitude spectrum
5. compares the numerical peaks with the theoretical normal-mode frequencies

This is a useful example because it links theory, simulation, and data analysis
in one compact workflow.

## Build-system tutorial

The `Makefile` is worth studying carefully because many research codes start
small and then grow. This project already shows a few good habits:

- source files are separated into `src/` and `include/`
- shared code is compiled once into a reusable object file
- build products and generated data go into dedicated directories
- plotting is automated instead of being an undocumented manual step

Targets:

- `make`
  Build the C programs.
- `make run`
  Build and execute the C programs, generating CSV data.
- `make plot`
  Regenerate data and build PNG figures with `gnuplot`.
- `make plot-python`
  Regenerate data and build PNG figures with Python `matplotlib`.
- `make clean`
  Remove build artifacts, generated data, plots, and local matplotlib cache
  files.

## Output files

Running `make run` creates:

- `output/good_sampling_signal.csv`
- `output/good_sampling_spectrum.csv`
- `output/undersampled_signal.csv`
- `output/undersampled_spectrum.csv`
- `output/short_record_signal.csv`
- `output/short_record_spectrum.csv`
- `output/coupled_oscillators_time.csv`
- `output/coupled_oscillators_spectrum.csv`

These CSV files are intentionally simple so that you can inspect them with any
tool you already know: Python, Julia, MATLAB, spreadsheets, or command-line
utilities.

## Plotting with gnuplot

Two ready-to-run `gnuplot` scripts are included:

- `scripts/plot_sampling.gp`
- `scripts/plot_coupled_oscillators.gp`

Run:

```bash
make plot
```

This:

1. rebuilds the programs if needed
2. regenerates the CSV data in `output/`
3. generates PNG figures in `plots/`

## Plotting with Python

Two Python scripts using `matplotlib` are included:

- `scripts/plot_sampling.py`
- `scripts/plot_coupled_oscillators.py`

Run:

```bash
make plot-python
```

This:

1. regenerates the CSV data in `output/`
2. runs the Python plotting scripts
3. writes PNG figures into `plots/`

The Python coupled-oscillator plot also overlays the theoretical normal-mode
frequencies as vertical reference lines. That is a good habit in research:
compare theory and numerics directly in the figure when possible.

## How to use this as a tutorial exercise

Suggested workflow:

1. Run `make run` and read the terminal output before looking at the plots.
2. Predict what the FFT should show for each case.
3. Open the CSV files and inspect their columns.
4. Generate plots and compare the visual evidence with the console output.
5. Change one parameter at a time and rerun the whole workflow.

If you want a more structured format, work through the guided worksheet in
`LESSONS.md`.

For instructors or for self-check after attempting the worksheet, concise model
answers are collected in `SOLUTIONS.md`.

## Suggested modifications for students

1. Change the sample rate in `src/sampling_demo.c` and predict the alias before
   rerunning the code.
2. Increase the duration of the short-record case until the 50 Hz and 55 Hz
   peaks separate clearly.
3. Remove the mean-subtraction step in `src/fft_utils.c` and inspect how the
   DC component changes the spectra.
4. Change the spring constants in `src/coupled_oscillators_fft.c` and compare
   the measured peaks with the updated theoretical values.
5. Change the initial conditions so that one normal mode is excited more
   strongly than the other.
6. Write a third physics example, such as a damped oscillator or a driven
   oscillator, and reuse the FFT helper code.

## Common beginner mistakes this project tries to prevent

- applying an FFT to non-uniformly sampled data
- confusing Nyquist frequency with frequency resolution
- assuming more samples always fix a too-short acquisition time
- plotting only the final figure without checking raw generated data
- treating the build and plotting steps as separate undocumented manual work

## Final remark

The project is deliberately modest in size, but the workflow is close to what
you will do in real research: define a model, generate data, document the
assumptions behind that data, analyze it with the right numerical tool, and
make the process reproducible.
