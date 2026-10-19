# Guided Worksheet

This short worksheet is meant for either:

- a supervised lab session
- self-study after reading the tutorial code

The goal is not only to run the programs, but to practice making predictions,
checking assumptions, and connecting numerical output to physical reasoning.

## Part 1: Before running anything

1. In your own words, what does an FFT tell you about a signal?
2. What is the difference between:
   - sampling rate
   - Nyquist frequency
   - frequency resolution
3. Why is a uniformly spaced time grid important for the FFT routines used in
   this project?

Write short answers before moving on.

## Part 2: Build and run

Run:

```bash
make run
```

Answer:

1. Which CSV files are created in `output/`?
2. Which program studies pure sampling issues?
3. Which program studies a mechanics problem?
4. What two frequencies are reported in the well-sampled case?

## Part 3: Sampling and aliasing

Focus on `src/sampling_demo.c`.

1. In the `undersampled` case, the sampling frequency is $f_s = 128\,\mathrm{Hz}$.
2. Using $f_N = \dfrac{f_s}{2}$, what is the Nyquist frequency $f_N$ for this case?
3. Why can a $120\,\mathrm{Hz}$ component not be reconstructed faithfully here?
4. The program reports an aliased peak near `8 Hz`.
5. Explain physically and mathematically why a high-frequency signal can appear
   as a false low-frequency signal after sampling.

Check your reasoning against the generated spectrum.

## Part 4: Frequency resolution

Focus on the `short_record` case.

1. The signal contains $50\,\mathrm{Hz}$ and $55\,\mathrm{Hz}$.
2. Using $T = N \Delta t$, what is the total record length $T$?
3. Using $\Delta f = \dfrac{f_s}{N} = \dfrac{1}{T}$, what is the frequency spacing $\Delta f$ of the FFT bins?
4. Why does this make the two frequencies difficult to separate?
5. If you wanted to improve the separation, would you first change:
   - the total acquisition time
   - the plotting tool
   - the output file format

Explain your choice.

## Part 5: Coupled oscillators

Focus on `src/coupled_oscillators_fft.c`.

1. Write down the equations of motion used in the code in the form
   $m \ddot{x}_1 = \cdots$ and $m \ddot{x}_2 = \cdots$.
2. What physical system do these equations represent?
3. What are the two normal modes of this system?
4. Why does the FFT of $x_1(t)$ show more than one important frequency?
5. Why is it useful to compare the numerical FFT peaks with the theoretical
   normal-mode frequencies?

## Part 6: Plot inspection

Generate plots with one of:

```bash
make plot
```

or

```bash
make plot-python
```

Answer:

1. Which plot makes aliasing easiest to see?
2. Which plot makes limited frequency resolution easiest to see?
3. In the coupled-oscillator FFT plot, do the numerical peaks line up well with
   the theoretical reference frequencies?
4. If they are not exactly identical, what numerical reasons could explain the
   small difference?

## Part 7: Small code modifications

Choose one of the following tasks.

### Option A: Change the sampling rate

Edit `src/sampling_demo.c` and modify the sampling frequency of one case.

Before rerunning:

1. Predict whether aliasing will occur.
2. Predict the approximate peak locations.

After rerunning:

3. Compare your prediction with the FFT output.

### Option B: Change the observation time

Edit the short-record case.

1. Increase the number of samples while keeping the interpretation clear.
2. Predict how $\Delta f = \dfrac{1}{T}$ will change.
3. Check whether the $50\,\mathrm{Hz}$ and $55\,\mathrm{Hz}$ peaks become more distinguishable.

### Option C: Change the spring constants

Edit `src/coupled_oscillators_fft.c`.

1. Change `wall_k` or `coupling_k`.
2. Predict how the two normal-mode frequencies
   $f_{\mathrm{in}} = \dfrac{1}{2\pi}\sqrt{\dfrac{k}{m}}$ and
   $f_{\mathrm{out}} = \dfrac{1}{2\pi}\sqrt{\dfrac{k + 2 k_c}{m}}$
   should change.
3. Rerun the program and compare the new FFT peaks with your prediction.

## Part 8: Reflection

Write a short paragraph addressing all three points:

1. One thing you learned about FFTs from this project.
2. One thing you learned about scientific software organization.
3. One thing you would improve or extend in the project.

## Optional extension

If you want a more research-like exercise, add one new example:

- a damped oscillator
- a driven oscillator
- a signal with added random noise

For your new example, state:

1. the model
2. the parameters
3. the expected spectral signature
4. whether the FFT confirms your expectation
