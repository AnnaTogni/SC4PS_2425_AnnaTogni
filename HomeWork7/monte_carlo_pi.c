// Exercise 2 -- Monte Carlo estimate of pi with two LCG streams
//
// Throw N random points (x, y) uniformly in the unit square [0,1) x [0,1).
// The fraction inside the quarter circle x^2 + y^2 <= 1 estimates its area pi/4:
//   pi_N = 4 * (points inside) / N.
// As in the notebook, x and y come from two LCG STREAMS: the same generator started
// from two different seeds (123 for x, 987654321 for y).
//
// Each "inside" is a Bernoulli variable with p = pi/4, so
//   sigma(pi_N) = 4 sqrt(p (1-p) / N) = sqrt(pi (4 - pi) / N) ~ 1.64 / sqrt(N).
// For every N the estimate is repeated R times (the streams simply continue) to
// measure the RMS error, which should follow this 1/sqrt(N) law.
//
// Usage: ./monte_carlo_pi [N_max] [R] [seed_x] [seed_y]
//        (default N_max = 1e7, R = 20, seeds 123 and 987654321)
// Output: output/monte_carlo_pi.csv

#include "common.h"

#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double estimate_pi(lcg_t *gx, lcg_t *gy, unsigned long N)
{
    unsigned long inside = 0;
    for (unsigned long i = 0; i < N; i++) {
        double x = lcg_uniform(gx);
        double y = lcg_uniform(gy);
        inside += (x * x + y * y <= 1.0);
    }
    return 4.0 * (double)inside / (double)N;
}

int main(int argc, char *argv[])
{
    unsigned long N_max = arg_ulong(argc, argv, 1, 10000000, "N_max");
    unsigned long R = arg_ulong(argc, argv, 2, 20, "R (repetitions)");
    uint32_t seed_x = arg_seed(argc, argv, 3, 123u);
    uint32_t seed_y = arg_seed(argc, argv, 4, 987654321u);
    lcg_t gx = lcg_new(seed_x); // stream 1: x coordinates
    lcg_t gy = lcg_new(seed_y); // stream 2: y coordinates

    FILE *f = open_output("monte_carlo_pi.csv");
    fprintf(f, "N,pi_first_run,error_first_run,rms_error,expected_sigma\n");

    printf("Exercise 2: Monte Carlo estimate of pi, two LCG streams (seeds %u, %u), "
           "%lu repetitions per N\n\n", seed_x, seed_y, R);
    printf("%11s %14s %12s %12s %14s %8s\n", "N", "pi (1st run)", "error", "RMS error",
           "sigma theory", "ratio");

    double sum_lx = 0, sum_ly = 0, sum_lxx = 0, sum_lxy = 0;
    int points = 0;
    for (unsigned long N = 10; N <= N_max; N *= 10) {
        double first = 0.0, sum_sq = 0.0;
        for (unsigned long k = 0; k < R; k++) {
            double est = estimate_pi(&gx, &gy, N);
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
    return 0;
}
