// exponential.h -- exponential sample by inverse transform (Exercises 4 and 5)
//
// If U ~ Uniform(0,1) and F(y) = 1 - exp(-lambda y) is the exponential CDF, then
//   Y = F^{-1}(U) = -ln(1 - U) / lambda
// has CDF F: P(Y <= y) = P(U <= F(y)) = F(y).
// gsl_rng_uniform returns U in [0, 1), so 1 - U is in (0, 1] and the log is always finite.
// (-ln(U)/lambda would work too, since 1 - U is also uniform, but ln(0) is possible.)

#ifndef EXPONENTIAL_H
#define EXPONENTIAL_H

#include <math.h>
#include <gsl/gsl_rng.h>

#define LAMBDA 1.5

static inline void exponential_sample(gsl_rng *r, double lambda, double *y, unsigned long n)
{
    for (unsigned long i = 0; i < n; i++)
        y[i] = -log(1.0 - gsl_rng_uniform(r)) / lambda;
}

#endif
