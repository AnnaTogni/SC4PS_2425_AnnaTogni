// Exercise 2 -- Monte Carlo estimate of pi with two RNG streams
//
// Throw N random points (x, y) uniformly in the unit square [0,1) x [0,1).
// The fraction inside the quarter circle x^2 + y^2 <= 1 estimates its area pi/4:
//   pi_N = 4 * (points inside) / N.
// As in the notebook, x and y come from two RNG STREAMS: the same generator (the LCG) started
// from two different seeds (123 for x, 987654321 for y).
//
// Each "inside" is a Bernoulli variable with p = pi/4, so the expected error is
//   sigma(pi_N) = 4 sqrt(p (1-p) / N) = sqrt(pi (4 - pi) / N) ~ 1.64 / sqrt(N).
// The program estimates pi for N = 10, 100, ..., N_max and compares the absolute error
// |pi_N - pi| with this 1/sqrt(N) law.
//
// Usage: ./monte_carlo_pi [N_max] [seed_x] [seed_y]
//        (default N_max = 1e7, seeds 123 and 987654321)
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
    uint32_t seed_x = arg_seed(argc, argv, 2, 123u);
    uint32_t seed_y = arg_seed(argc, argv, 3, 987654321u);
    lcg_t gx = lcg_new(seed_x); // stream 1: x coordinates
    lcg_t gy = lcg_new(seed_y); // stream 2: y coordinates

    FILE *f = open_output("monte_carlo_pi.csv");
    fprintf(f, "N,pi_estimate,abs_error,expected_sigma\n");

    printf("Exercise 2: Monte Carlo estimate of pi, two RNG streams (seeds %u, %u)\n\n", seed_x,
           seed_y);
    printf("%11s %14s %12s %14s\n", "N", "pi estimate", "|error|", "sigma theory");

    for (unsigned long N = 10; N <= N_max; N *= 10) {
        double est = estimate_pi(&gx, &gy, N);
        double err = fabs(est - M_PI);
        double sigma = sqrt(M_PI * (4.0 - M_PI) / (double)N);
        printf("%11lu %14.8f %12.2e %14.2e\n", N, est, err, sigma);
        fprintf(f, "%lu,%.12f,%.6e,%.6e\n", N, est, err, sigma);
        fflush(stdout);
    }
    fclose(f);

    printf("\nResults written to %s/monte_carlo_pi.csv\n", OUTPUT_DIR);
    return 0;
}
