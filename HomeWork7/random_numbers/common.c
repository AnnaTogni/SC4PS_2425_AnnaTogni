// common.c -- helpers shared by the HomeWork7 programs (see common.h)

#define _POSIX_C_SOURCE 200809L // mkdir() with -std=c11

#include "common.h"

#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>

#include <gsl/gsl_cdf.h>

gsl_rng *rng_new(unsigned long seed)
{
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    if (r == NULL) {
        fprintf(stderr, "Error: cannot allocate the random number generator\n");
        exit(1);
    }
    gsl_rng_set(r, seed);
    return r;
}

unsigned long arg_ulong(int argc, char *argv[], int i, unsigned long def, const char *name)
{
    if (i >= argc)
        return def;
    char *end;
    double v = strtod(argv[i], &end);
    if (end == argv[i] || *end != '\0' || v < 1 || v != floor(v) || v > 4e18) {
        fprintf(stderr, "Error: %s must be a positive integer, got '%s'\n", name, argv[i]);
        exit(1);
    }
    return (unsigned long)v;
}

FILE *open_output(const char *name)
{
    if (mkdir(OUTPUT_DIR, 0775) != 0 && errno != EEXIST) {
        perror(OUTPUT_DIR);
        exit(1);
    }
    char path[256];
    snprintf(path, sizeof(path), "%s/%s", OUTPUT_DIR, name);
    FILE *f = fopen(path, "w");
    if (f == NULL) {
        perror(path);
        exit(1);
    }
    return f;
}

double chi_square(const size_t *observed, const double *expected, size_t bins, size_t *dof,
                  double *p_value)
{
    double chi2 = 0.0, obs = 0.0, exp = 0.0;
    size_t merged = 0;
    for (size_t k = 0; k < bins; k++) {
        obs += (double)observed[k];
        exp += expected[k];
        // close the merged bin when it has enough expected counts (or at the end)
        if (exp >= 5.0 || k == bins - 1) {
            if (exp > 0.0) {
                chi2 += (obs - exp) * (obs - exp) / exp;
                merged++;
            }
            obs = exp = 0.0;
        }
    }
    *dof = merged > 1 ? merged - 1 : 1;
    *p_value = gsl_cdf_chisq_Q(chi2, (double)*dof);
    return chi2;
}
