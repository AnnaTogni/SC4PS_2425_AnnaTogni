// common.h -- helpers shared by the HomeWork7 programs
//
// Random numbers: GSL generator "mt19937" (Mersenne Twister), period 2^19937 - 1.
// A fixed seed makes every run reproducible; a different seed gives an independent sample.

#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stddef.h>

#include <gsl/gsl_rng.h>

#define OUTPUT_DIR   "output"
#define DEFAULT_SEED 12345UL

// Allocate a Mersenne Twister generator with the given seed (exits on failure).
gsl_rng *rng_new(unsigned long seed);

// Read an optional positive integer argument argv[i] (e.g. a sample size or a seed);
// accepts 1e6. Returns def if i >= argc. Exits with a message on invalid input.
unsigned long arg_ulong(int argc, char *argv[], int i, unsigned long def, const char *name);

// Create OUTPUT_DIR if needed and open OUTPUT_DIR/name for writing (exits on failure).
FILE *open_output(const char *name);

// Chi-square test of a histogram against expected counts.
// Bins with expected count < 5 are merged with the next one (standard rule).
// Returns chi2 and sets *dof (degrees of freedom = merged bins - 1) and *p_value.
double chi_square(const size_t *observed, const double *expected, size_t bins, size_t *dof,
                  double *p_value);

#endif
