// Exercise 3 -- Change of variables: Y = U^2 with U ~ Uniform(0,1)
//
// For y in (0,1):  F_Y(y) = P(U^2 <= y) = P(U <= sqrt(y)) = sqrt(y)
//                  f_Y(y) = dF_Y/dy = 1 / (2 sqrt(y))      (diverges at y = 0, but integrable)
//
// The histogram of Y is compared with f_Y. Since f_Y changes a lot inside the first bins,
// the fair comparison is with the AVERAGE density in each bin,
//   (F_Y(b) - F_Y(a)) / (b - a) = (sqrt(b) - sqrt(a)) / (b - a),
// not with f_Y at the bin centre. A chi-square test compares observed and expected counts.
//
// Usage: ./square_uniform [n] [bins] [seed]   (default n = 1e6, bins = 50, seed = 12345)
// Output: output/square_uniform.csv

#include "common.h"

#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    unsigned long n = arg_ulong(argc, argv, 1, 1000000, "n");
    size_t bins = arg_ulong(argc, argv, 2, 50, "bins");
    unsigned long seed = arg_ulong(argc, argv, 3, DEFAULT_SEED, "seed");
    gsl_rng *r = rng_new(seed);

    size_t *count = calloc(bins, sizeof(*count));
    double *expected = malloc(bins * sizeof(*expected));
    if (count == NULL || expected == NULL) {
        fprintf(stderr, "Error: not enough memory\n");
        return 1;
    }

    double sum = 0.0, sum2 = 0.0;
    for (unsigned long i = 0; i < n; i++) {
        double u = gsl_rng_uniform(r); // [0, 1)
        double y = u * u;              // the change of variables
        sum += y;
        sum2 += y * y;
        size_t k = (size_t)(y * bins);
        count[k < bins ? k : bins - 1]++;
    }

    FILE *f = open_output("square_uniform.csv");
    fprintf(f, "bin_left,bin_right,bin_center,count,density_histogram,density_bin_average,"
               "pdf_at_center\n");
    double width = 1.0 / (double)bins;
    for (size_t k = 0; k < bins; k++) {
        double a = k * width, b = (k + 1) * width, c = 0.5 * (a + b);
        double p_bin = sqrt(b) - sqrt(a); // F(b) - F(a)
        expected[k] = n * p_bin;
        fprintf(f, "%.6f,%.6f,%.6f,%zu,%.8f,%.8f,%.8f\n", a, b, c, count[k],
                count[k] / (n * width), p_bin / width, 1.0 / (2.0 * sqrt(c)));
    }
    fclose(f);

    size_t dof;
    double p_value;
    double chi2 = chi_square(count, expected, bins, &dof, &p_value);
    double mean = sum / n, var = sum2 / n - mean * mean;

    printf("Exercise 3: Y = U^2, n = %lu samples, %zu bins (seed %lu)\n\n", n, bins, seed);
    printf("                 sample        exact\n");
    printf("mean     E[Y]  = %.6f   1/3  = %.6f\n", mean, 1.0 / 3.0);
    printf("variance Var[Y]= %.6f   4/45 = %.6f\n", var, 4.0 / 45.0);
    printf("first bin [0, %.3f): histogram density %.4f, average of f_Y %.4f, f_Y(center) %.4f\n",
           width, count[0] / (n * width), sqrt(width) / width, 1.0 / (2.0 * sqrt(0.5 * width)));
    printf("\nchi-square = %.1f with %zu degrees of freedom (expected ~ %zu +- %.0f), p-value = %.3f\n",
           chi2, dof, dof, sqrt(2.0 * dof), p_value);
    printf("-> %s\n", p_value > 0.01 ? "histogram consistent with f_Y(y) = 1/(2 sqrt(y))"
                                     : "histogram NOT consistent with f_Y (p < 0.01)");
    printf("Histogram written to %s/square_uniform.csv\n", OUTPUT_DIR);

    free(count);
    free(expected);
    gsl_rng_free(r);
    return 0;
}
