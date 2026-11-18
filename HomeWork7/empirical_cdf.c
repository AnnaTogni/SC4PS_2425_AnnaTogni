// Exercise 5 -- Empirical CDF of the exponential sample of Exercise 4
//
// With the same n and seed as ./exponential, exponential_sample() produces exactly the
// same sample. Sorted, y_(1) <= ... <= y_(n), the empirical CDF is the step function
//   F_n(y) = (number of samples <= y) / n,   which jumps by 1/n at every y_(i).
// It is compared with the exact CDF F(y) = 1 - exp(-lambda y).
//
// Usage: ./empirical_cdf [n] [seed]     (default n = 200000, seed = 246813579, as Exercise 4)
// Output: output/empirical_cdf.csv (at most ~2000 points of the step function, for the plot)

#include "common.h"
#include "exponential.h"

#include <stdlib.h>

#define PLOT_POINTS 2000

static double exact_cdf(double y) { return 1.0 - exp(-LAMBDA * y); }

static int compare_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
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
    qsort(y, n, sizeof(*y), compare_double);

    printf("Exercise 5: empirical CDF of the exponential sample (lambda = %.2f, n = %lu, seed %u)\n\n",
           LAMBDA, n, seed);
    printf("%10s %12s %12s\n", "y", "F_n(y)", "F(y)");
    const double show[] = {0.1, 0.5, 1.0, 2.0, 3.0};
    for (size_t k = 0; k < sizeof(show) / sizeof(show[0]); k++) {
        unsigned long below = 0; // number of samples <= show[k] (y is sorted)
        while (below < n && y[below] <= show[k])
            below++;
        printf("%10.2f %12.6f %12.6f\n", show[k], (double)below / n, exact_cdf(show[k]));
    }

    // Points of the step function for the plot: every k-th order statistic
    FILE *f = open_output("empirical_cdf.csv");
    fprintf(f, "y,ecdf,exact_cdf\n");
    unsigned long step = n > PLOT_POINTS ? n / PLOT_POINTS : 1;
    for (unsigned long i = 0; i < n; i += step)
        fprintf(f, "%.8f,%.8f,%.8f\n", y[i], (double)(i + 1) / n, exact_cdf(y[i]));
    fprintf(f, "%.8f,%.8f,%.8f\n", y[n - 1], 1.0, exact_cdf(y[n - 1]));
    fclose(f);
    printf("\nEmpirical CDF written to %s/empirical_cdf.csv\n", OUTPUT_DIR);

    free(y);
    return 0;
}
