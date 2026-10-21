// Bonus -- RANDU's planes (the notebook's warning example) vs the notebook LCG
//
// RANDU: x_{n+1} = 65539 x_n mod 2^31, u_n = x_n / 2^31.
// Since 65539 = 2^16 + 3, squaring gives 65539^2 = 6 * 65539 - 9 (mod 2^31), so
//     x_{n+2} = 6 x_{n+1} - 9 x_n   (mod 2^31)
// and therefore, dividing by 2^31,
//     9 u_n - 6 u_{n+1} + u_{n+2} = k,   k an INTEGER.
// Every triple (u_n, u_{n+1}, u_{n+2}) lies on one of the parallel planes
// 9x - 6y + z = k; since the left side is between -6 and 10, only 15 planes are hit.
//
// For a good generator, s = 9 u_n - 6 u_{n+1} + u_{n+2} is a continuous variable: the
// fraction of triples with s within 0.01 of an integer should be about 2%.
//
// Usage: ./randu_planes [n] [seed]     (default n = 1e6 triples, seed = 1)
// Output: output/randu_triples.csv, output/lcg_triples.csv (3000 triples each, for the plot)
//         output/plane_combination.csv (s for 100000 triples of each generator)

#include "common.h"

#include <math.h>
#include <stdlib.h>

#define PLOT_TRIPLES 3000
#define HIST_TRIPLES 100000

int main(int argc, char *argv[])
{
    unsigned long n = arg_ulong(argc, argv, 1, 1000000, "n");
    uint32_t seed = arg_seed(argc, argv, 2, 1u);

    uint32_t xr = seed | 1u; // RANDU needs an odd seed
    lcg_t g = lcg_new(seed);
    double r0 = randu_next(&xr) / 2147483648.0, r1 = randu_next(&xr) / 2147483648.0;
    double l0 = lcg_uniform(&g), l1 = lcg_uniform(&g);

    FILE *fr = open_output("randu_triples.csv");
    FILE *fl = open_output("lcg_triples.csv");
    FILE *fs = open_output("plane_combination.csv");
    fprintf(fr, "u0,u1,u2\n");
    fprintf(fl, "u0,u1,u2\n");
    fprintf(fs, "s_randu,s_lcg\n");

    unsigned long near_r = 0, near_l = 0;
    long kmin = 100, kmax = -100;
    double worst = 0.0;
    for (unsigned long i = 0; i < n; i++) {
        double r2 = randu_next(&xr) / 2147483648.0;
        double l2 = lcg_uniform(&g);
        double sr = 9.0 * r0 - 6.0 * r1 + r2;
        double sl = 9.0 * l0 - 6.0 * l1 + l2;

        long k = lround(sr);
        worst = fmax(worst, fabs(sr - k)); // 0 up to rounding if RANDU's planes are exact
        kmin = k < kmin ? k : kmin;
        kmax = k > kmax ? k : kmax;
        near_r += fabs(sr - k) < 0.01;
        near_l += fabs(sl - lround(sl)) < 0.01;

        if (i < PLOT_TRIPLES) {
            fprintf(fr, "%.10f,%.10f,%.10f\n", r0, r1, r2);
            fprintf(fl, "%.10f,%.10f,%.10f\n", l0, l1, l2);
        }
        if (i < HIST_TRIPLES)
            fprintf(fs, "%.10f,%.10f\n", sr, sl);
        r0 = r1;
        r1 = r2;
        l0 = l1;
        l1 = l2;
    }
    fclose(fr);
    fclose(fl);
    fclose(fs);

    printf("Bonus: correlations of consecutive triples (u_n, u_{n+1}, u_{n+2}), %lu triples\n\n", n);
    printf("s = 9 u_n - 6 u_{n+1} + u_{n+2}\n");
    printf("RANDU: s is always an integer (largest distance from an integer: %.1e);\n"
           "       it takes %ld values (%ld ... %ld): all triples lie on %ld planes\n",
           worst, kmax - kmin + 1, kmin, kmax, kmax - kmin + 1);
    printf("fraction of triples with s within 0.01 of an integer (uniform triples: ~2%%):\n");
    printf("  RANDU          %6.2f%%\n", 100.0 * near_r / n);
    printf("  notebook LCG   %6.2f%%\n", 100.0 * near_l / n);
    printf("Triples written to %s/randu_triples.csv, %s/lcg_triples.csv\n", OUTPUT_DIR, OUTPUT_DIR);
    return 0;
}
