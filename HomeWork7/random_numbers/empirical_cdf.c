// Exercise 5 -- Empirical CDF of the exponential sample of Exercise 4
//
// With the same n and seed as ./exponential, exponential_sample() produces exactly the
// same sample. Sorted, y_(1) <= ... <= y_(n), the empirical CDF is the step function
//   F_n(y) = (number of samples <= y) / n,   which jumps by 1/n at every y_(i).
// It is compared with the exact CDF F(y) = 1 - exp(-lambda y) through the
// Kolmogorov-Smirnov distance
//   D_n = max_y |F_n(y) - F(y)| = max_i max( i/n - F(y_(i)), F(y_(i)) - (i-1)/n ).
// For a correct sample, sqrt(n) D_n follows the Kolmogorov distribution: it is below
// 1.36 in 95% of the cases, so D_n shrinks like 1/sqrt(n).
//
// Usage: ./empirical_cdf [n] [seed]     (default n = 200000, seed = 246813579, as Exercise 4)
// Output: output/empirical_cdf.csv (at most ~2000 points of the step function, for the plot)
//         output/ks_vs_n.csv       (D_n for n = 10, 100, ..., n)

#include "common.h"
#include "exponential.h"

#include <stdlib.h>

#include <gsl/gsl_sort.h>

#define PLOT_POINTS 2000

static double exact_cdf(double y) { return 1.0 - exp(-LAMBDA * y); }

// Kolmogorov-Smirnov distance of a SORTED sample from the exact CDF.
static double ks_distance(const double *y, unsigned long n, double *where)
{
    double d = 0.0;
    for (unsigned long i = 0; i < n; i++) {
        double F = exact_cdf(y[i]);
        double above = (double)(i + 1) / n - F; // F_n just after the jump
        double below = F - (double)i / n;       // F_n just before the jump
        double m = above > below ? above : below;
        if (m > d) {
            d = m;
            *where = y[i];
        }
    }
    return d;
}

int main(int argc, char *argv[])
{
    unsigned long n = arg_ulong(argc, argv, 1, 200000, "n");
    uint32_t seed = arg_seed(argc, argv, 2, EXPONENTIAL_SEED);

    double *y = malloc(n * sizeof(*y));
    if (y == NULL) {
        fprintf(stderr, "Error: not enough memory\n");
        return 1;
    }
    lcg_t gen = lcg_new(seed);
    exponential_sample(&gen, LAMBDA, y, n); // same sample as Exercise 4

    // D_n for the first 10, 100, ... values of the same sample: how fast F_n -> F
    FILE *g = open_output("ks_vs_n.csv");
    fprintf(g, "n,D_n,sqrt_n_D_n\n");
    printf("Exercise 5: empirical CDF of the exponential sample (lambda = %.2f, seed %u)\n\n",
           LAMBDA, seed);
    printf("%10s %12s %12s\n", "n", "D_n", "sqrt(n) D_n");
    double *tmp = malloc(n * sizeof(*tmp));
    for (unsigned long m = 10; tmp != NULL && m <= n; m *= 10) {
        for (unsigned long i = 0; i < m; i++)
            tmp[i] = y[i];
        gsl_sort(tmp, 1, m);
        double where;
        double d = ks_distance(tmp, m, &where);
        printf("%10lu %12.6f %12.3f\n", m, d, sqrt((double)m) * d);
        fprintf(g, "%lu,%.8f,%.6f\n", m, d, sqrt((double)m) * d);
    }
    free(tmp);
    fclose(g);

    gsl_sort(y, 1, n);
    double where = 0.0;
    double d = ks_distance(y, n, &where);
    printf("\nFull sample: n = %lu, D_n = %.2e at y = %.3f, sqrt(n) D_n = %.3f\n", n, d, where,
           sqrt((double)n) * d);
    printf("95%% critical value: sqrt(n) D_n < 1.36  -> %s\n",
           sqrt((double)n) * d < 1.36 ? "the sample is consistent with the exact CDF"
                                      : "the sample is NOT consistent with the exact CDF");

    // Points of the step function for the plot: every k-th order statistic
    FILE *f = open_output("empirical_cdf.csv");
    fprintf(f, "y,ecdf,exact_cdf\n");
    unsigned long step = n > PLOT_POINTS ? n / PLOT_POINTS : 1;
    for (unsigned long i = 0; i < n; i += step)
        fprintf(f, "%.8f,%.8f,%.8f\n", y[i], (double)(i + 1) / n, exact_cdf(y[i]));
    fprintf(f, "%.8f,%.8f,%.8f\n", y[n - 1], 1.0, exact_cdf(y[n - 1]));
    fclose(f);
    printf("Written %s/empirical_cdf.csv and %s/ks_vs_n.csv\n", OUTPUT_DIR, OUTPUT_DIR);

    free(y);
    return 0;
}
