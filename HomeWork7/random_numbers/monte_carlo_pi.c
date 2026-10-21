// Exercise 2 -- Monte Carlo estimate of pi with two RNG streams
//
// Throw N random points (x, y) uniformly in the unit square [0,1) x [0,1).
// The fraction inside the quarter circle x^2 + y^2 < 1 estimates its area pi/4:
//   pi_N = 4 * (points inside) / N.
// x and y come from two INDEPENDENT generator streams (two Mersenne Twisters with
// different seeds), so the two coordinates of a point are not correlated.
//
// Each "inside" is a Bernoulli variable with p = pi/4, so
//   sigma(pi_N) = 4 sqrt(p (1-p) / N) = sqrt(pi (4 - pi) / N) ~ 1.64 / sqrt(N).
// For every N the estimate is repeated R times (new random numbers each time) to
// measure the RMS error, which should follow this 1/sqrt(N) law.
//
// Usage: ./monte_carlo_pi [N_max] [R] [seed]   (default N_max = 1e7, R = 20, seed = 12345)
// Output: output/monte_carlo_pi.csv

#include "common.h"

#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double estimate_pi(gsl_rng *rx, gsl_rng *ry, unsigned long N)
{
    unsigned long inside = 0;
    for (unsigned long i = 0; i < N; i++) {
        double x = gsl_rng_uniform(rx);
        double y = gsl_rng_uniform(ry);
        inside += (x * x + y * y < 1.0);
    }
    return 4.0 * (double)inside / (double)N;
}

int main(int argc, char *argv[])
{
    unsigned long N_max = arg_ulong(argc, argv, 1, 10000000, "N_max");
    unsigned long R = arg_ulong(argc, argv, 2, 20, "R (repetitions)");
    unsigned long seed = arg_ulong(argc, argv, 3, DEFAULT_SEED, "seed");
    gsl_rng *rx = rng_new(seed);     // stream 1: x coordinates
    gsl_rng *ry = rng_new(seed + 1); // stream 2: y coordinates

    FILE *f = open_output("monte_carlo_pi.csv");
    fprintf(f, "N,pi_first_run,error_first_run,rms_error,expected_sigma\n");

    printf("Exercise 2: Monte Carlo estimate of pi, %lu repetitions per N (seeds %lu, %lu)\n\n",
           R, seed, seed + 1);
    printf("%11s %14s %12s %12s %14s %8s\n", "N", "pi (1st run)", "error", "RMS error",
           "sigma theory", "ratio");

    double sum_lx = 0, sum_ly = 0, sum_lxx = 0, sum_lxy = 0;
    int points = 0;
    for (unsigned long N = 10; N <= N_max; N *= 10) {
        double first = 0.0, sum_sq = 0.0;
        for (unsigned long k = 0; k < R; k++) {
            double est = estimate_pi(rx, ry, N);
            if (k == 0)
                first = est;
            sum_sq += (est - M_PI) * (est - M_PI);
        }
        double rms = sqrt(sum_sq / (double)R);
        double sigma = sqrt(M_PI * (4.0 - M_PI) / (double)N);
        printf("%11lu %14.8f %12.2e %12.2e %14.2e %8.2f\n", N, first, first - M_PI, rms, sigma,
               rms / sigma);
        fprintf(f, "%lu,%.12f,%.6e,%.6e,%.6e\n", N, first, first - M_PI, rms, sigma);
        fflush(stdout);

        // least-squares fit of log10(rms) = a + b log10(N)
        double lx = log10((double)N), ly = log10(rms);
        sum_lx += lx;
        sum_ly += ly;
        sum_lxx += lx * lx;
        sum_lxy += lx * ly;
        points++;
    }
    fclose(f);

    double slope = (points * sum_lxy - sum_lx * sum_ly) / (points * sum_lxx - sum_lx * sum_lx);
    printf("\nFit RMS error ~ N^b: b = %.3f   (theory: -0.5, i.e. 100x more points -> 10x smaller"
           " error)\n", slope);
    printf("Results written to %s/monte_carlo_pi.csv\n", OUTPUT_DIR);
    gsl_rng_free(rx);
    gsl_rng_free(ry);
    return 0;
}
