// Exercise 1 -- Coin tosses and the Law of Large Numbers
//
// N fair coin tosses with the notebook LCG: heads if u < 0.5 (as in the notebook).
// The running fraction of heads
//   f_n = (number of heads in the first n tosses) / n
// converges to 1/2 (law of large numbers). Its typical distance from 1/2 is the
// standard deviation of f_n: sigma_n = sqrt(p (1-p) / n) = 1 / (2 sqrt(n)).
//
// u < 0.5 looks at the HIGHEST bit of the LCG state. For comparison the program also
// "tosses" with the LOWEST bit (x mod 2): in an LCG mod 2^32 that bit simply
// alternates 0, 1, 0, 1, ... -- a fraction that is "too good" to be random.
//
// Usage: ./coin_tosses [N] [seed]        (default N = 100000, seed = 123456789)
// Output: output/coin_tosses.csv  (n, running fraction, sigma_n, low-bit fraction)

#include "common.h"

#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    unsigned long N = arg_ulong(argc, argv, 1, 100000, "N");
    uint32_t seed = arg_seed(argc, argv, 2, 123456789u);
    lcg_t g = lcg_new(seed);
    FILE *f = open_output("coin_tosses.csv");
    fprintf(f, "n,fraction_heads,sigma,fraction_low_bit\n");

    printf("Exercise 1: %lu fair coin tosses with the notebook LCG (seed %u)\n\n", N, seed);
    printf("%10s %10s %12s %12s %10s   %s\n", "n", "heads", "fraction", "|f - 0.5|", "in sigmas",
           "low-bit 'coin'");

    unsigned long heads = 0, low_heads = 0, next_print = 10;
    for (unsigned long n = 1; n <= N; n++) {
        uint32_t x = lcg_next(&g);
        heads += (x / 4294967296.0 < 0.5); // u < 0.5: the highest bit of x is 0
        low_heads += x & 1u;               // the lowest bit of x
        double frac = (double)heads / (double)n;
        double sigma = 0.5 / sqrt((double)n);
        fprintf(f, "%lu,%.10f,%.10f,%.10f\n", n, frac, sigma, (double)low_heads / (double)n);
        if (n == next_print || n == N) {
            printf("%10lu %10lu %12.6f %12.6f %10.2f   %.6f\n", n, heads, frac, fabs(frac - 0.5),
                   fabs(frac - 0.5) / sigma, (double)low_heads / (double)n);
            if (n == next_print)
                next_print *= 10;
        }
    }
    fclose(f);

    lcg_t h = lcg_new(seed);
    printf("\nLowest bit of the first 16 states: ");
    for (int i = 0; i < 16; i++)
        printf("%u", lcg_next(&h) & 1u);
    printf("  <- strictly alternating: never use the low bits of an LCG!\n");
    printf("sigma_n = 1/(2 sqrt(n)) is the expected size of |f_n - 0.5|: the error shrinks only\n"
           "like 1/sqrt(n) (10x more tosses -> ~3.2x smaller error).\n");
    printf("Running fraction written to %s/coin_tosses.csv\n", OUTPUT_DIR);
    return 0;
}
