// HomeWork5 -- Floating-point summation
//
// Sums a vector of doubles in several ways and compares the results:
//
//   1. plain for loop
//   2. GSL gsl_vector_sum()
//   3. Kahan compensated summation
//   4. Neumaier compensated summation ("improved Kahan")
//   5. pairwise (recursive halving) summation
//   6. for loop after sorting by increasing |value|
//   7. for loop in long double (more mantissa bits)
//   8. exact sum, correctly rounded to double (Shewchuk's algorithm, as Python's math.fsum)
//
// Default vector (lecture exercise): {1.0, 1.0e16, -1.0e16, -0.5}, exact sum = 0.5.
//
// ******************************************************
// HOW TO COMPILE:
//   make
// HOW TO RUN:
//   ./floating_point_sum                 (default vector)
//   ./floating_point_sum 1e16 1 -1e16    (any list of numbers)
// ******************************************************

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_vector.h>

#define MAX_TRACE 10 // print the step-by-step trace only for short vectors

// ****************************************************************************
// Summation methods
// ****************************************************************************

static double sum_loop(const double *v, size_t n)
{
    double s = 0.0;
    for (size_t i = 0; i < n; i++)
        s += v[i];
    return s;
}

static double sum_gsl(double *v, size_t n)
{
    gsl_vector_view view = gsl_vector_view_array(v, n);
    return gsl_vector_sum(&view.vector);
}

// Kahan: c keeps the low-order part lost in each addition.
// The formula c = (t - s) - y recovers the lost part exactly ONLY if |s| >= |y|.
static double sum_kahan(const double *v, size_t n, int trace)
{
    double s = 0.0, c = 0.0;
    if (trace)
        printf("  Kahan    | %-10s | %-22s | %-22s\n", "add", "running sum s", "compensation c");
    for (size_t i = 0; i < n; i++) {
        double y = v[i] - c;
        double t = s + y;
        c = (t - s) - y;
        s = t;
        if (trace)
            printf("  step %zu   | %-10g | %-22.17g | %-22.17g\n", i + 1, v[i], s, c);
    }
    return s;
}

// Neumaier: like Kahan, but computes the lost part from the LARGER of the two
// operands, so it is exact also when the new term is bigger than the sum.
// The correction is added only once, at the end.
static double sum_neumaier(const double *v, size_t n, int trace)
{
    double s = 0.0, c = 0.0;
    if (trace)
        printf("  Neumaier | %-10s | %-22s | %-22s\n", "add", "running sum s", "correction c");
    for (size_t i = 0; i < n; i++) {
        double t = s + v[i];
        if (fabs(s) >= fabs(v[i]))
            c += (s - t) + v[i];
        else
            c += (v[i] - t) + s;
        s = t;
        if (trace)
            printf("  step %zu   | %-10g | %-22.17g | %-22.17g\n", i + 1, v[i], s, c);
    }
    return s + c;
}

// Pairwise: sum the two halves recursively, then add them.
static double sum_pairwise(const double *v, size_t n)
{
    if (n == 0)
        return 0.0;
    if (n == 1)
        return v[0];
    size_t half = n / 2;
    return sum_pairwise(v, half) + sum_pairwise(v + half, n - half);
}

static int compare_abs(const void *a, const void *b)
{
    double x = fabs(*(const double *)a), y = fabs(*(const double *)b);
    return (x > y) - (x < y);
}

static double sum_sorted(const double *v, size_t n)
{
    double *w = malloc(n * sizeof(*w));
    if (w == NULL)
        return NAN;
    for (size_t i = 0; i < n; i++)
        w[i] = v[i];
    qsort(w, n, sizeof(*w), compare_abs);
    double s = sum_loop(w, n);
    free(w);
    return s;
}

static double sum_long_double(const double *v, size_t n)
{
    long double s = 0.0L;
    for (size_t i = 0; i < n; i++)
        s += v[i];
    return (double)s;
}

// Exact sum, correctly rounded (Shewchuk 1997, the algorithm of Python's math.fsum).
// The running sum is kept EXACTLY as a list of non-overlapping doubles ("partials");
// each new value is added with error-free transformations (hi + lo = a + b exactly).
static double sum_exact(const double *v, size_t n)
{
    double *p = malloc((n + 1) * sizeof(*p));
    size_t np = 0;
    if (p == NULL)
        return NAN;

    for (size_t i = 0; i < n; i++) {
        double x = v[i];
        size_t k = 0;
        for (size_t j = 0; j < np; j++) {
            double y = p[j];
            if (fabs(x) < fabs(y)) {
                double tmp = x;
                x = y;
                y = tmp;
            }
            double hi = x + y;
            double lo = y - (hi - x); // exact rounding error of x + y
            if (lo != 0.0)
                p[k++] = lo;
            x = hi;
        }
        p[k++] = x;
        np = k;
    }

    // Add the partials from the largest; stop when the rest can no longer matter.
    double hi = 0.0;
    if (np > 0) {
        size_t j = np - 1;
        hi = p[j];
        while (j > 0) {
            double x = hi, y = p[--j];
            hi = x + y;
            double lo = y - (hi - x);
            if (lo != 0.0) {
                // Round-half-even correction, as in math.fsum
                if (j > 0 && ((lo < 0 && p[j - 1] < 0) || (lo > 0 && p[j - 1] > 0))) {
                    double y2 = lo * 2.0, x2 = hi + y2;
                    if (y2 == x2 - hi)
                        hi = x2;
                }
                break;
            }
        }
    }
    free(p);
    return hi;
}

// ****************************************************************************

int main(int argc, char *argv[])
{
    static const double default_v[] = {1.0, 1.0e16, -1.0e16, -0.5};
    size_t n = argc > 1 ? (size_t)(argc - 1) : sizeof(default_v) / sizeof(default_v[0]);
    double *v = malloc(n * sizeof(*v));
    if (v == NULL) {
        fprintf(stderr, "Error: not enough memory\n");
        return 1;
    }

    for (size_t i = 0; i < n; i++) {
        if (argc > 1) {
            char *end;
            v[i] = strtod(argv[i + 1], &end);
            if (end == argv[i + 1] || *end != '\0') {
                fprintf(stderr, "Error: '%s' is not a number\n"
                                "Usage: %s [x1 x2 ...]   (default: 1 1e16 -1e16 -0.5)\n",
                        argv[i + 1], argv[0]);
                free(v);
                return 1;
            }
        } else {
            v[i] = default_v[i];
        }
    }

    printf("Vector (%zu values):", n);
    for (size_t i = 0; i < n && i < 20; i++)
        printf(" %g", v[i]);
    printf("%s\n\n", n > 20 ? " ..." : "");

    int trace = n <= MAX_TRACE;
    if (trace) {
        printf("Step by step:\n");
        sum_kahan(v, n, 1);
        printf("\n");
        sum_neumaier(v, n, 1);
        printf("\n");
    }

    double exact = sum_exact(v, n);
    struct {
        const char *name;
        double value;
    } results[] = {
        {"1. for loop", sum_loop(v, n)},
        {"2. GSL gsl_vector_sum", sum_gsl(v, n)},
        {"3. Kahan", sum_kahan(v, n, 0)},
        {"4. Neumaier", sum_neumaier(v, n, 0)},
        {"5. pairwise", sum_pairwise(v, n)},
        {"6. sorted by |x|, then loop", sum_sorted(v, n)},
        {"7. long double loop", sum_long_double(v, n)},
        {"8. exact (Shewchuk / fsum)", exact},
    };

    printf("%-30s %-24s %s\n", "Method", "Result", "Correct?");
    printf("%-30s %-24s %s\n", "------", "------", "--------");
    for (size_t i = 0; i < sizeof(results) / sizeof(results[0]); i++)
        printf("%-30s %-24.17g %s\n", results[i].name, results[i].value,
               results[i].value == exact ? "YES" : "NO");

    free(v);
    return 0;
}
