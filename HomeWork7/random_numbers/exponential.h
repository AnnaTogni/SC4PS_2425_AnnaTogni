// exponential.h -- exponential sample by inverse transform (Exercises 4 and 5)
//
// If U ~ Uniform(0,1) and F(y) = 1 - exp(-lambda y) is the exponential CDF, then
//   Y = F^{-1}(U) = -ln(1 - U) / lambda
// has CDF F: P(Y <= y) = P(U <= F(y)) = F(y).
// The LCG gives U = x / 2^32 in [0, 1), so 1 - U is in (0, 1] and the log is always finite.
// (-ln(U)/lambda would work too, since 1 - U is also uniform, but U = 0 would give ln(0)).
// Largest possible value: U = 1 - 2^-32 -> Y = 32 ln(2) / lambda ~ 14.8.

#ifndef EXPONENTIAL_H
#define EXPONENTIAL_H

#include <math.h>

#include "common.h"

#define LAMBDA 1.5

#define EXPONENTIAL_SEED 246813579u // seed of the notebook's Exercise 4 solution

static inline void exponential_sample(lcg_t *g, double lambda, double *y, unsigned long n)
{
    for (unsigned long i = 0; i < n; i++)
        y[i] = -log(1.0 - lcg_uniform(g)) / lambda;
}

#endif
