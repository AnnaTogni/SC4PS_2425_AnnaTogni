/*
 * coupled_oscillators_fft.c
 *
 * Tutorial purpose:
 * Connect a standard mechanics problem to a numerical-analysis workflow.
 *
 * The program:
 *
 * 1. defines the equations of motion for two coupled oscillators
 * 2. integrates them with a GSL ODE solver
 * 3. samples the displacement on a uniform time grid
 * 4. computes an FFT of the sampled displacement
 * 5. compares the numerical peaks with the theoretical normal-mode frequencies
 *
 * This example is useful for physics students because it shows that a numerical
 * pipeline often has several stages, each with its own assumptions.
 */

#include "fft_utils.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    double mass;
    double wall_k;
    double coupling_k;
} SystemParameters;

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

static int rhs(double t, const double y[], double f[], void *params)
{
    const SystemParameters *p = params;
    const double x1 = y[0];
    const double v1 = y[1];
    const double x2 = y[2];
    const double v2 = y[3];

    (void)t;

    /*
     * State vector convention:
     *
     * y[0] = x1
     * y[1] = v1
     * y[2] = x2
     * y[3] = v2
     *
     * Two identical masses attached to walls and coupled by a spring:
     *
     * m x1'' = -(k + kc) x1 + kc x2
     * m x2'' =  kc x1 - (k + kc) x2
     */
    f[0] = v1;
    f[1] = (-(p->wall_k + p->coupling_k) * x1 + p->coupling_k * x2) / p->mass;
    f[2] = v2;
    f[3] = (p->coupling_k * x1 - (p->wall_k + p->coupling_k) * x2) / p->mass;

    return GSL_SUCCESS;
}

static int write_time_series_csv(const char *path,
                                 const double *times,
                                 const double *x1,
                                 const double *x2,
                                 size_t n)
{
    FILE *fp = fopen(path, "w");
    size_t i;

    if (fp == NULL) {
        return -1;
    }

    fprintf(fp, "time_s,x1_m,x2_m\n");
    for (i = 0; i < n; ++i) {
        fprintf(fp, "%.12f,%.12f,%.12f\n", times[i], x1[i], x2[i]);
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

static void report_mode_near_frequency(const double *frequencies,
                                       const double *amplitudes,
                                       size_t bins,
                                       double target_frequency_hz,
                                       const char *label)
{
    size_t start;
    size_t stop;
    size_t best_bin;
    size_t i;

    /*
     * Search a few bins around the expected theoretical value. The FFT bins are
     * discrete, so the measured peak is generally near the theoretical value,
     * not exactly equal to it.
     */
    best_bin = (size_t)llround(target_frequency_hz / (frequencies[1] - frequencies[0]));
    start = best_bin > 3U ? best_bin - 3U : 1U;
    stop = best_bin + 3U < bins ? best_bin + 3U : bins - 1U;

    for (i = start; i <= stop; ++i) {
        if (amplitudes[i] > amplitudes[best_bin]) {
            best_bin = i;
        }
    }

    printf("  %s: expected %.5f Hz, measured %.5f Hz, amplitude %.4f\n",
           label,
           target_frequency_hz,
           frequencies[best_bin],
           amplitudes[best_bin]);
}

int main(int argc, char **argv)
{
    const char *output_dir = argc > 1 ? argv[1] : "output";
    const SystemParameters params = {
        .mass = 1.0,
        .wall_k = 25.0,
        .coupling_k = 7.0,
    };
    const size_t sample_count = 1024U;
    const double dt = 0.125;
    const size_t bins = sample_count / 2U + 1U;
    const double omega_in_phase = sqrt(params.wall_k / params.mass);
    const double omega_out_of_phase =
        sqrt((params.wall_k + 2.0 * params.coupling_k) / params.mass);
    const double f_in_phase = omega_in_phase / (2.0 * M_PI);
    const double f_out_of_phase = omega_out_of_phase / (2.0 * M_PI);
    double *times = NULL;
    double *x1 = NULL;
    double *x2 = NULL;
    double *frequencies = NULL;
    double *amplitudes = NULL;
    char time_series_path[256];
    char spectrum_path[256];
    gsl_odeiv2_system system = {rhs, NULL, 4, (void *)&params};
    gsl_odeiv2_driver *driver = NULL;
    double t = 0.0;
    double y[4] = {0.10, 0.0, 0.0, 0.0};
    size_t i;
    int status = EXIT_FAILURE;

    if (!fftu_is_power_of_two(sample_count)) {
        fprintf(stderr, "sample_count must be a power of two\n");
        return EXIT_FAILURE;
    }
    if (ensure_directory(output_dir) != 0) {
        fprintf(stderr, "cannot create output directory '%s'\n", output_dir);
        return EXIT_FAILURE;
    }

    times = malloc(sample_count * sizeof(*times));
    x1 = malloc(sample_count * sizeof(*x1));
    x2 = malloc(sample_count * sizeof(*x2));
    frequencies = malloc(bins * sizeof(*frequencies));
    amplitudes = malloc(bins * sizeof(*amplitudes));
    if (times == NULL || x1 == NULL || x2 == NULL ||
        frequencies == NULL || amplitudes == NULL) {
        fprintf(stderr, "memory allocation failed\n");
        goto cleanup;
    }

    /*
     * The integrator uses adaptive internal time steps for accuracy, but we
     * request the solution on a uniform output grid. That combination is useful
     * in practice: ODE solvers stay accurate and the FFT still receives equally
     * spaced samples.
     */
    driver = gsl_odeiv2_driver_alloc_y_new(&system,
                                           gsl_odeiv2_step_rkf45,
                                           1e-6,
                                           1e-9,
                                           1e-9);
    if (driver == NULL) {
        fprintf(stderr, "failed to allocate GSL ODE driver\n");
        goto cleanup;
    }

    times[0] = 0.0;
    x1[0] = y[0];
    x2[0] = y[2];

    for (i = 1; i < sample_count; ++i) {
        const double target_time = i * dt;
        const int gsl_status = gsl_odeiv2_driver_apply(driver, &t, target_time, y);

        if (gsl_status != GSL_SUCCESS) {
            fprintf(stderr,
                    "ODE integration failed at sample %zu with status %d\n",
                    i,
                    gsl_status);
            goto cleanup;
        }

        /* Store the solution only on the uniform grid used for the FFT. */
        times[i] = target_time;
        x1[i] = y[0];
        x2[i] = y[2];
    }

    if (fftu_compute_one_sided_amplitude_spectrum(x1,
                                                  sample_count,
                                                  dt,
                                                  frequencies,
                                                  amplitudes) != 0) {
        fprintf(stderr, "FFT analysis failed\n");
        goto cleanup;
    }

    snprintf(time_series_path,
             sizeof(time_series_path),
             "%s/coupled_oscillators_time.csv",
             output_dir);
    snprintf(spectrum_path,
             sizeof(spectrum_path),
             "%s/coupled_oscillators_spectrum.csv",
             output_dir);

    if (write_time_series_csv(time_series_path, times, x1, x2, sample_count) != 0) {
        fprintf(stderr, "failed to write %s\n", time_series_path);
        goto cleanup;
    }
    if (write_spectrum_csv(spectrum_path, frequencies, amplitudes, bins) != 0) {
        fprintf(stderr, "failed to write %s\n", spectrum_path);
        goto cleanup;
    }

    printf("Coupled oscillators analyzed with an FFT\n");
    printf("=======================================\n");
    printf("mass = %.2f kg, wall spring = %.2f N/m, coupling spring = %.2f N/m\n",
           params.mass,
           params.wall_k,
           params.coupling_k);
    printf("sample spacing dt = %.3f s, total time = %.3f s, df = %.5f Hz\n",
           dt,
           sample_count * dt,
           1.0 / (sample_count * dt));
    printf("theoretical formulas: f_in = sqrt(k/m)/(2 pi), f_out = sqrt((k + 2 kc)/m)/(2 pi)\n");
    printf("normal-mode frequencies from theory:\n");
    report_mode_near_frequency(frequencies, amplitudes, bins, f_in_phase, "in-phase mode");
    report_mode_near_frequency(frequencies, amplitudes, bins, f_out_of_phase, "out-of-phase mode");
    printf("wrote %s\n", time_series_path);
    printf("wrote %s\n", spectrum_path);

    status = EXIT_SUCCESS;

cleanup:
    gsl_odeiv2_driver_free(driver);
    free(times);
    free(x1);
    free(x2);
    free(frequencies);
    free(amplitudes);
    return status;
}
