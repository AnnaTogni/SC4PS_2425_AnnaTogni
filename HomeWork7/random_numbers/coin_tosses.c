// Exercise 1 -- Coin tosses and the Law of Large Numbers
//
// N fair coin tosses with the notebook LCG: heads if u < 0.5 (as in the notebook).
// The running fraction of heads
//   f_n = (number of heads in the first n tosses) / n
// converges to 1/2 (law of large numbers). Its typical distance from 1/2 is the
// standard deviation of f_n: sigma_n = sqrt(p (1-p) / n) = 1 / (2 sqrt(n)).
//
// Usage: ./coin_tosses [N] [seed]        (default N = 100000, seed = 123456789)
// Output: output/coin_tosses.csv  (n, running fraction, sigma_n)

#include "common.h"

#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    unsigned long N = arg_ulong(argc, argv, 1, 100000, "N");
    uint32_t seed = arg_seed(argc, argv, 2, 123456789u);
    lcg_t g = lcg_new(seed);
    FILE *f = open_output("coin_tosses.csv");
    fprintf(f, "n,fraction_heads,sigma\n");

    printf("Exercise 1: %lu fair coin tosses with the notebook LCG (seed %u)\n\n", N, seed);
    printf("%10s %10s %12s %12s %10s\n", "n", "heads", "fraction", "|f - 0.5|", "in sigmas");

    unsigned long heads = 0, next_print = 10;
    for (unsigned long n = 1; n <= N; n++) {
        uint32_t x = lcg_next(&g);
        heads += (x / 4294967296.0 < 0.5); // heads if u < 0.5
        double frac = (double)heads / (double)n;
        double sigma = 0.5 / sqrt((double)n);
        fprintf(f, "%lu,%.10f,%.10f\n", n, frac, sigma);
        if (n == next_print || n == N) {
            printf("%10lu %10lu %12.6f %12.6f %10.2f\n", n, heads, frac, fabs(frac - 0.5),
                   fabs(frac - 0.5) / sigma);
            if (n == next_print)
                next_print *= 10;
        }
    }
    fclose(f);

    printf("\nsigma_n = 1/(2 sqrt(n)) is the expected size of |f_n - 0.5|: the error shrinks only\n"
           "like 1/sqrt(n) (10x more tosses -> ~3.2x smaller error).\n");
    printf("Running fraction written to %s/coin_tosses.csv\n", OUTPUT_DIR);
    return 0;
}
