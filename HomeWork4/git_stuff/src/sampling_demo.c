/*
 * sampling_demo.c
 *
 * Tutorial purpose:
 * Demonstrate how sampling choices affect what the FFT can recover from a
 * signal. This is one of the most important beginner lessons in numerical data
 * analysis for physics students.
 *
 * The program studies three cases:
 *
 * 1. good_sampling
 *    The sampling rate is high enough to capture both frequency components.
 *
 * 2. undersampled
 *    The 120 Hz component lies above the Nyquist frequency and appears as an
 *    alias at a false lower frequency.
 *
 * 3. short_record
 *    The record length is too short to separate 50 Hz and 55 Hz clearly,
 *    illustrating limited frequency resolution.
 *
 * The generated CSV files can be plotted with either gnuplot or Python.
 */

#include "fft_utils.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef double (*signal_fn)(double t);

/*
 * Each SamplingCase defines one self-contained numerical experiment.
 *
 * name            file prefix for generated CSV output
 * lesson          short explanation printed to the terminal
 * signal          function to sample
 * sample_rate_hz  samples per second
 * sample_count    number of uniformly spaced samples
 */
typedef struct {
    const char *name;
    const char *lesson;
    signal_fn signal;
    double sample_rate_hz;
    size_t sample_count;
} SamplingCase;

static double aliasing_signal(double t)
{
    /* Two clean sinusoidal components used to reveal aliasing clearly. */
    return sin(2.0 * M_PI * 50.0 * t) + 0.70 * sin(2.0 * M_PI * 120.0 * t);
}

static double close_frequencies_signal(double t)
{
    /* Two nearby frequencies used to illustrate limited frequency resolution. */
    return sin(2.0 * M_PI * 50.0 * t) + 0.85 * sin(2.0 * M_PI * 55.0 * t);
}

static int ensure_directory(const char *path)
{
    struct stat st;

    if (stat(path, &st) == 0) {
        return S_ISDIR(st.st_mode) ? 0 : -1;
    }

    if (mkdir(path, 0775) == 0) {
        return 0;
    }

    return errno == EEXIST ? 0 : -1;
}

static int write_signal_csv(const char *path,
                            const double *times,
                            const double *samples,
                            size_t n)
{
    FILE *fp = fopen(path, "w");
    size_t i;

    if (fp == NULL) {
        return -1;
    }

    fprintf(fp, "time_s,signal\n");
    for (i = 0; i < n; ++i) {
        fprintf(fp, "%.12f,%.12f\n", times[i], samples[i]);
    }

    fclose(fp);
    return 0;
}

static int write_spectrum_csv(const char *path,
                              const double *frequencies,
                              const double *amplitudes,
                              size_t bins)
{
    FILE *fp = fopen(path, "w");
    size_t i;

    if (fp == NULL) {
        return -1;
    }

    fprintf(fp, "frequency_hz,amplitude\n");
    for (i = 0; i < bins; ++i) {
        fprintf(fp, "%.12f,%.12f\n", frequencies[i], amplitudes[i]);
    }

    fclose(fp);
    return 0;
}

static void insert_peak(size_t bin,
                        double amplitude,
                        size_t *peak_bins,
                        double *peak_amplitudes,
                        size_t peak_count)
{
    size_t i;

    /* Maintain a tiny sorted list of the strongest candidate peaks. */
    for (i = 0; i < peak_count; ++i) {
        if (amplitude > peak_amplitudes[i]) {
            size_t j;

            for (j = peak_count - 1; j > i; --j) {
                peak_amplitudes[j] = peak_amplitudes[j - 1];
                peak_bins[j] = peak_bins[j - 1];
            }
            peak_amplitudes[i] = amplitude;
            peak_bins[i] = bin;
            break;
        }
    }
}

static void report_dominant_peaks(const double *frequencies,
                                  const double *amplitudes,
                                  size_t bins)
{
    size_t peak_bins[3] = {0, 0, 0};
    double peak_amplitudes[3] = {-1.0, -1.0, -1.0};
    double max_amplitude = 0.0;
    const double relative_threshold = 1e-3;
    size_t i;

    for (i = 1; i < bins; ++i) {
        if (amplitudes[i] > max_amplitude) {
            max_amplitude = amplitudes[i];
        }
    }

    for (i = 1; i + 1 < bins; ++i) {
        if (amplitudes[i] >= amplitudes[i - 1] &&
            amplitudes[i] > amplitudes[i + 1]) {
            /*
             * Identify simple local maxima. For a teaching example we prefer a
             * clear and readable rule instead of a complicated peak finder.
             */
            insert_peak(i, amplitudes[i], peak_bins, peak_amplitudes, 3);
        }
    }

    for (i = 0; i < 3; ++i) {
        if (peak_amplitudes[i] > relative_threshold * max_amplitude) {
            printf("    peak %zu: f = %8.3f Hz, amplitude = %.4f\n",
                   i + 1,
                   frequencies[peak_bins[i]],
                   peak_amplitudes[i]);
        }
    }
}

static int analyze_case(const SamplingCase *sample_case, const char *output_dir)
{
    const double dt = 1.0 / sample_case->sample_rate_hz;
    const size_t bins = sample_case->sample_count / 2U + 1U;
    const double nyquist_hz = 0.5 * sample_case->sample_rate_hz;
    const double resolution_hz =
        sample_case->sample_rate_hz / (double)sample_case->sample_count;
    double *times = NULL;
    double *samples = NULL;
    double *frequencies = NULL;
    double *amplitudes = NULL;
    char signal_path[256];
    char spectrum_path[256];
    size_t i;
    int status = 0;

    times = malloc(sample_case->sample_count * sizeof(*times));
    samples = malloc(sample_case->sample_count * sizeof(*samples));
    frequencies = malloc(bins * sizeof(*frequencies));
    amplitudes = malloc(bins * sizeof(*amplitudes));
    if (times == NULL || samples == NULL || frequencies == NULL ||
        amplitudes == NULL) {
        status = -1;
        goto cleanup;
    }

    for (i = 0; i < sample_case->sample_count; ++i) {
        /* Uniform sampling: t_i = i dt. */
        times[i] = i * dt;
        samples[i] = sample_case->signal(times[i]);
    }

    status = fftu_compute_one_sided_amplitude_spectrum(samples,
                                                       sample_case->sample_count,
                                                       dt,
                                                       frequencies,
                                                       amplitudes);
    if (status != 0) {
        goto cleanup;
    }

    snprintf(signal_path,
             sizeof(signal_path),
             "%s/%s_signal.csv",
             output_dir,
             sample_case->name);
    snprintf(spectrum_path,
             sizeof(spectrum_path),
             "%s/%s_spectrum.csv",
             output_dir,
             sample_case->name);

    if (write_signal_csv(signal_path, times, samples, sample_case->sample_count) !=
        0) {
        status = -2;
        goto cleanup;
    }
    if (write_spectrum_csv(spectrum_path, frequencies, amplitudes, bins) != 0) {
        status = -3;
        goto cleanup;
    }

    printf("%s\n", sample_case->name);
    printf("  %s\n", sample_case->lesson);
    printf("  sampling rate = %.1f Hz, Nyquist = %.1f Hz, duration = %.3f s\n",
           sample_case->sample_rate_hz,
           nyquist_hz,
           sample_case->sample_count * dt);
    printf("  FFT frequency resolution = %.3f Hz\n", resolution_hz);
    printf("  dominant spectral peaks:\n");
    report_dominant_peaks(frequencies, amplitudes, bins);
    printf("  wrote %s\n", signal_path);
    printf("  wrote %s\n\n", spectrum_path);

cleanup:
    free(times);
    free(samples);
    free(frequencies);
    free(amplitudes);
    return status;
}

int main(int argc, char **argv)
{
    const char *output_dir = argc > 1 ? argv[1] : "output";
    const SamplingCase cases[] = {
        {
            "good_sampling",
            "The 50 Hz and 120 Hz components are both below Nyquist, so the FFT recovers both frequencies correctly.",
            aliasing_signal,
            512.0,
            512U,
        },
        {
            "undersampled",
            "The 120 Hz component is above Nyquist for fs = 128 Hz and folds back to an alias near 8 Hz.",
            aliasing_signal,
            128.0,
            128U,
        },
        {
            "short_record",
            "The signal contains 50 Hz and 55 Hz, but the record is too short, so the FFT cannot separate nearby frequencies cleanly.",
            close_frequencies_signal,
            512.0,
            64U,
        },
    };
    size_t i;

    if (ensure_directory(output_dir) != 0) {
        fprintf(stderr, "cannot create output directory '%s'\n", output_dir);
        return EXIT_FAILURE;
    }

    if (!fftu_is_power_of_two(cases[0].sample_count)) {
        fprintf(stderr, "internal error: sample counts must be powers of two\n");
        return EXIT_FAILURE;
    }

    printf("Sampling and FFT demonstration\n");
    printf("==============================\n");
    printf("This program writes CSV files that can be plotted with gnuplot, Python, or Julia.\n");
    printf("Read the three cases as a small numerical experiment on sampling strategy.\n\n");

    for (i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
        if (!fftu_is_power_of_two(cases[i].sample_count)) {
            fprintf(stderr, "case '%s' does not use a power-of-two sample count\n", cases[i].name);
            return EXIT_FAILURE;
        }
        if (analyze_case(&cases[i], output_dir) != 0) {
            fprintf(stderr, "failed to analyze case '%s'\n", cases[i].name);
            return EXIT_FAILURE;
        }
    }

    printf("Interpretation\n");
    printf("--------------\n");
    printf("1. A sampling rate above twice the highest frequency avoids aliasing.\n");
    printf("2. A longer time window improves frequency resolution because df = fs / N = 1 / T.\n");
    printf("3. FFTs assume evenly spaced samples; non-uniform data require different tools.\n");

    return EXIT_SUCCESS;
}
