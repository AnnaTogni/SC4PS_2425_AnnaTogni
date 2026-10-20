// HomeWork2
//
// Addition of two vectors : d = a*x + y
//   a     : scalar
//   x, y  : vectors of dimension N with all elements equal
//
// Everything is read from the command line, the result is checked
// element by element, and an optional expected value (e.g. 7.4) is
// compared both exactly (==) and with a floating-point tolerance.
//
// ******************************************************
// HOW TO COMPILE:
//   make                      (or: gcc -Wall -Wextra -O2 NumericalSum.c -o numericalsum)
//
// HOW TO RUN:
//   ./numericalsum N a x y [expected]
//   ./numericalsum 10  3 0.1 7.1 7.4
//   ./numericalsum 1e6 3 0.1 7.1 7.4
//   ./numericalsum 1e8 3 0.1 7.1 7.4
// ******************************************************

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Relative tolerance for the "approximately equal" test: a few units of
// machine precision (DBL_EPSILON ~ 2.2e-16 is the gap between 1.0 and the
// next representable double).
#define REL_TOL (4.0 * DBL_EPSILON)

// Number of elements printed at the beginning of d (the rest is only checked).
#define N_PRINT 5

static void usage(const char *prog)
{
    fprintf(stderr,
            "Usage: %s N a x y [expected]\n"
            "  N         vector dimension (e.g. 10, 1e6, 1e8)\n"
            "  a         scalar\n"
            "  x, y      value of every element of x and y\n"
            "  expected  optional value that every d[i] should be equal to\n"
            "Example: %s 1e6 3 0.1 7.1 7.4\n",
            prog, prog);
}

// Read a double from a string; return 0 if the string is not a valid number.
static int parse_double(const char *s, double *out)
{
    char *end;
    *out = strtod(s, &end);
    return end != s && *end == '\0';
}

// Return 1 if a and b are equal within a relative tolerance.
static int almost_equal(double a, double b)
{
    return fabs(a - b) <= REL_TOL * fmax(fabs(a), fabs(b));
}

static double seconds_now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

int main(int argc, char *argv[])
{
    if (argc != 5 && argc != 6) {
        usage(argv[0]);
        return 1;
    }

    // ---- 1. Read the input from the command line --------------------------
    // N is read as a double so that "1e8" is accepted, then checked to be a
    // positive integer.
    double N_in, a, xv, yv, expected = 0.0;
    int have_expected = (argc == 6);

    if (!parse_double(argv[1], &N_in) || N_in < 1 || N_in != floor(N_in)) {
        fprintf(stderr, "Error: N must be a positive integer, got '%s'\n", argv[1]);
        return 1;
    }
    if (!parse_double(argv[2], &a) || !parse_double(argv[3], &xv) ||
        !parse_double(argv[4], &yv) ||
        (have_expected && !parse_double(argv[5], &expected))) {
        fprintf(stderr, "Error: a, x, y and expected must be numbers\n");
        usage(argv[0]);
        return 1;
    }
    size_t N = (size_t)N_in;

    printf("N = %zu, a = %g, x = %g, y = %g\n", N, a, xv, yv);
    printf("Memory needed: 3 vectors x %zu x %zu bytes = %.1f MB\n",
           N, sizeof(double), 3.0 * N * sizeof(double) / 1e6);
    fflush(stdout);   // show the lines above before a possible error message

    // ---- 2. Allocate the vectors (and check that it worked) ---------------
    double *x = malloc(N * sizeof(double));
    double *y = malloc(N * sizeof(double));
    double *d = malloc(N * sizeof(double));
    if (x == NULL || y == NULL || d == NULL) {
        fprintf(stderr, "Error: not enough memory for N = %zu\n", N);
        free(x);
        free(y);
        free(d);
        return 1;
    }

    for (size_t i = 0; i < N; i++) {
        x[i] = xv;
        y[i] = yv;
    }

    // ---- 3. Compute d = a*x + y ------------------------------------------
    double t0 = seconds_now();
    for (size_t i = 0; i < N; i++)
        d[i] = a * x[i] + y[i];
    double t1 = seconds_now();

    printf("Time for d = a*x + y: %.6f s\n\n", t1 - t0);

    printf("First elements of d:\n");
    for (size_t i = 0; i < N && i < N_PRINT; i++)
        printf("  d[%zu] = %.17g\n", i, d[i]);
    if (N > N_PRINT)
        printf("  ... (%zu more)\n", N - N_PRINT);

    // ---- 4. Check every element against the value computed by the program -
    double computed = a * xv + yv;
    size_t n_wrong = 0;
    for (size_t i = 0; i < N; i++)
        if (d[i] != computed)
            n_wrong++;

    printf("\nCheck 1: every d[i] == a*x + y = %.17g ?  %s (%zu of %zu differ)\n",
           computed, n_wrong == 0 ? "YES" : "NO", n_wrong, N);

    // ---- 5. Compare with the value the user expects (e.g. 7.4) ------------
    int ok = (n_wrong == 0);
    if (have_expected) {
        size_t n_exact = 0, n_close = 0;
        for (size_t i = 0; i < N; i++) {
            if (d[i] == expected)
                n_exact++;
            if (almost_equal(d[i], expected))
                n_close++;
        }
        printf("\nCheck 2: expected value %g (stored as %.17g)\n", expected, expected);
        printf("  exact comparison   d[i] == %g            : %s (%zu of %zu)\n",
               expected, n_exact == N ? "YES" : "NO", n_exact, N);
        printf("  with tolerance     |d[i] - %g| <= %.1e*|%g| : %s (%zu of %zu)\n",
               expected, REL_TOL, expected, n_close == N ? "YES" : "NO", n_close, N);
        printf("  difference d - expected = %.3e\n", computed - expected);
        if (n_exact != N && n_close == N)
            printf("  -> equal only within rounding: %g cannot be represented exactly in binary\n",
                   expected);
        ok = ok && (n_close == N);
    }

    free(x);
    free(y);
    free(d);

    return ok ? 0 : 2;
}
