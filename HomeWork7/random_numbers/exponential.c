// Exercise 4 -- Exponential distribution by inverse transform
//
//   Y = -ln(1 - U) / lambda,  U ~ Uniform(0,1),  lambda = 1.5
//   pdf  f(y) = lambda exp(-lambda y),   mean 1/lambda,   variance 1/lambda^2
//
// The histogram is compared with the pdf (bin-average density and chi-square test).
//
// U comes from the notebook LCG (see common.h).
//
// Usage: ./exponential [n] [bins] [seed]   (default n = 200000, bins = 60, seed = 246813579)
// Output: output/exponential.csv

#include "common.h"
#include "exponential.h"

#include <stdlib.h>

#define Y_MAX 6.0 // histogram range [0, Y_MAX); P(Y >= 6) = exp(-9) ~ 1e-4

int main(int argc, char *argv[])
{
    unsigned long n = arg_ulong(argc, argv, 1, 200000, "n");
    size_t bins = arg_ulong(argc, argv, 2, 60, "bins");
    uint32_t seed = arg_seed(argc, argv, 3, EXPONENTIAL_SEED);
    const double lambda = LAMBDA;

    double *y = malloc(n * sizeof(*y));
    size_t *count = calloc(bins + 1, sizeof(*count)); // last entry: overflow (y >= Y_MAX)
    double *expected = malloc((bins + 1) * sizeof(*expected));
    if (y == NULL || count == NULL || expected == NULL) {
        fprintf(stderr, "Error: not enough memory\n");
        return 1;
    }

    lcg_t g = lcg_new(seed);
    exponential_sample(&g, lambda, y, n);

    double width = Y_MAX / (double)bins, sum = 0.0, sum2 = 0.0, ymax = 0.0;
    for (unsigned long i = 0; i < n; i++) {
        sum += y[i];
        sum2 += y[i] * y[i];
        ymax = y[i] > ymax ? y[i] : ymax;
        size_t k = (size_t)(y[i] / width);
        count[k < bins ? k : bins]++;
    }

    FILE *f = open_output("exponential.csv");
    fprintf(f, "bin_left,bin_right,bin_center,count,density_histogram,density_bin_average,"
               "pdf_at_center\n");
    for (size_t k = 0; k < bins; k++) {
        double a = k * width, b = (k + 1) * width, c = 0.5 * (a + b);
        double p_bin = exp(-lambda * a) - exp(-lambda * b); // F(b) - F(a)
        expected[k] = n * p_bin;
        fprintf(f, "%.6f,%.6f,%.6f,%zu,%.8f,%.8f,%.8f\n", a, b, c, count[k],
                count[k] / (n * width), p_bin / width, lambda * exp(-lambda * c));
    }
    expected[bins] = n * exp(-lambda * Y_MAX); // overflow bin, included in the chi-square
    fclose(f);

    size_t dof;
    double p_value;
    double chi2 = chi_square(count, expected, bins + 1, &dof, &p_value);
    double mean = sum / n, var = sum2 / n - mean * mean;

    printf("Exercise 4: exponential, lambda = %.2f, n = %lu samples (notebook LCG, seed %u)\n\n",
           lambda, n, seed);
    printf("                  sample      exact\n");
    printf("mean     1/lambda   = %.5f   %.5f\n", mean, 1.0 / lambda);
    printf("variance 1/lambda^2 = %.5f   %.5f\n", var, 1.0 / (lambda * lambda));
    printf("largest value       = %.3f    (values >= %.0f: %zu, expected %.1f)\n", ymax, Y_MAX,
           count[bins], expected[bins]);
    printf("\nchi-square = %.1f with %zu degrees of freedom (expected ~ %zu +- %.0f), p-value = %.3f\n",
           chi2, dof, dof, sqrt(2.0 * dof), p_value);
    printf("-> %s\n", p_value > 0.01 ? "histogram consistent with f(y) = lambda exp(-lambda y)"
                                     : "histogram NOT consistent with the pdf (p < 0.01)");
    printf("Histogram written to %s/exponential.csv\n", OUTPUT_DIR);

    free(y);
    free(count);
    free(expected);
    return 0;
}
