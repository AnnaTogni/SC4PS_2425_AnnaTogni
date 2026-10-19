#ifndef FFT_UTILS_H
#define FFT_UTILS_H

#include <stddef.h>

/*
 * Public interface for the FFT helper routines used by this tutorial project.
 *
 * Why is this file separate from the application programs?
 *
 * In scientific software, the same numerical operation often appears in
 * multiple applications. Placing the declarations here and the implementation
 * in src/fft_utils.c keeps the code modular and reusable.
 *
 * The routines below assume uniformly sampled data. That assumption is not a
 * minor implementation detail; it is one of the central lessons of the
 * tutorial. The radix-2 FFT used here is appropriate when:
 *
 * - the number of samples is a power of two
 * - the sampling interval dt is constant
 * - the user wants the frequency content of a discrete time series
 */

/* Return nonzero if n is a power of two, zero otherwise. */
int fftu_is_power_of_two(size_t n);

/*
 * Compute the one-sided amplitude spectrum of a real signal.
 *
 * samples      input array of length n
 * n            number of samples, must be a power of two
 * dt           constant sampling step in seconds
 * frequencies  output array of length n / 2 + 1
 * amplitudes   output array of length n / 2 + 1
 *
 * Practical note for students:
 *
 * A real-valued time series has a symmetric complex Fourier transform. For that
 * reason, one often reports only the non-negative frequencies, known here as
 * the one-sided spectrum.
 *
 * The DC component and Nyquist bin use a scale factor of 1. All other bins are
 * multiplied by 2 so that a sine wave with amplitude A produces a spectral
 * peak close to A when it lands cleanly on a Fourier bin.
 *
 * Returns 0 on success and a negative value on error.
 */
int fftu_compute_one_sided_amplitude_spectrum(const double *samples,
                                              size_t n,
                                              double dt,
                                              double *frequencies,
                                              double *amplitudes);

#endif
