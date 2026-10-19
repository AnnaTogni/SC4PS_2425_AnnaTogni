/*
 * Shared FFT helper functions for the tutorial applications.
 *
 * This file is intentionally simple. Its job is to hide the low-level GSL FFT
 * call behind an interface that is easier to reason about in the teaching
 * examples. The main scientific ideas to notice are:
 *
 * - the sample count must match the requirements of the chosen FFT algorithm
 * - the input is treated as uniformly sampled data
 * - the mean is subtracted before the transform to reduce the DC component
 * - the output is converted into a one-sided amplitude spectrum
 */

#include "fft_utils.h"

#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_fft_complex.h>

int fftu_is_power_of_two(size_t n)
{
    return n > 0 && (n & (n - 1U)) == 0;
}

int fftu_compute_one_sided_amplitude_spectrum(const double *samples,
                                              size_t n,
                                              double dt,
                                              double *frequencies,
                                              double *amplitudes)
{
    double *data = NULL;
    double mean = 0.0;
    size_t i = 0;

    if (samples == NULL || frequencies == NULL || amplitudes == NULL) {
        return -1;
    }
    if (n < 2 || !fftu_is_power_of_two(n) || dt <= 0.0) {
        return -2;
    }

    data = calloc(2U * n, sizeof(*data));
    if (data == NULL) {
        return -3;
    }

    for (i = 0; i < n; ++i) {
        mean += samples[i];
    }
    mean /= (double)n;

    /*
     * Remove the mean before the FFT so that the DC component does not
     * dominate the spectrum when the interesting physics is in the oscillatory
     * part of the signal.
     */
    for (i = 0; i < n; ++i) {
        data[2U * i] = samples[i] - mean;
        data[2U * i + 1U] = 0.0;
    }

    if (gsl_fft_complex_radix2_forward(data, 1, n) != 0) {
        free(data);
        return -4;
    }

    /*
     * Convert the complex FFT result into a one-sided amplitude spectrum.
     *
     * The bin frequency is i / (n dt). For real input data, positive and
     * negative frequencies carry mirrored information, so interior bins are
     * multiplied by 2 to preserve the physical amplitude.
     */
    for (i = 0; i <= n / 2U; ++i) {
        const double re = data[2U * i];
        const double im = data[2U * i + 1U];
        const double scale = (i == 0 || i == n / 2U) ? 1.0 : 2.0;

        frequencies[i] = (double)i / ((double)n * dt);
        amplitudes[i] = scale * hypot(re, im) / (double)n;
    }

    free(data);
    return 0;
}
