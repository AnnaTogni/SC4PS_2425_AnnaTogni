// HomeWork5 -- Stability of Legendre polynomials
//
// For order m = 0 and degrees l = 0 ... l_max, compute the Legendre polynomials P_l(x)
// at several x (|x| < 1) in three ways:
//
//   1. forward three-term recurrence (Bonnet), in double precision:
//          P_0 = 1, P_1 = x,
//          P_{l+1} = [ (2l+1) x P_l - l P_{l-1} ] / (l+1)
//   2. high-precision reference: the same recurrence in quad precision
//          (__float128, 113-bit mantissa, about 34 decimal digits)
//   3. Miller-style backward recurrence, in double precision:
//          P~_{L+1} = 0, P~_L = 1,
//          P~_{l-1} = [ (2l+1) x P~_l - (l+1) P~_{l+1} ] / l,   l = L, L-1, ..., 1,
//          then P_l^(back) = c P~_l with c = 1 / P~_0
//
// and compares 1 and 3 with the reference (absolute and relative errors).
// GSL (gsl_sf_legendre_Pl_array) is also printed as a cross-check.
//
// ******************************************************
// HOW TO COMPILE:
//   make
// HOW TO RUN:
//   ./legendre_stability [l_max] [L] [x1 x2 ...]
//   ./legendre_stability                  (l_max = 50, L = 80, x = 0.1 0.5 0.9 0.99)
//   make plot                             (figures in plots/)
// ******************************************************

#define _POSIX_C_SOURCE 200809L // mkdir() with -std=c11

#include <errno.h>
#include <float.h>
#include <math.h>
#include <quadmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>

#define OUTPUT_DIR    "output"
#define DEFAULT_L_MAX 50
#define DEFAULT_L     80

static const double DEFAULT_X[] = {0.1, 0.5, 0.9, 0.99};

typedef __float128 quad;

// ****************************************************************************
// Recurrences
// ****************************************************************************

// 1. Forward recurrence from P_0 = 1, P_1 = x up to P_lmax.
static void forward_double(double x, int lmax, double *p)
{
    p[0] = 1.0;
    p[1] = x;
    for (int l = 1; l < lmax; l++)
        p[l + 1] = ((2.0 * l + 1.0) * x * p[l] - l * p[l - 1]) / (l + 1.0);
}

// 2. The same forward recurrence in quad precision: the reference values.
static void forward_quad(quad x, int lmax, quad *p)
{
    p[0] = 1;
    p[1] = x;
    for (int l = 1; l < lmax; l++)
        p[l + 1] = ((quad)(2 * l + 1) * x * p[l] - (quad)l * p[l - 1]) / (quad)(l + 1);
}

// 3. Miller-style backward recurrence: start from P~_{L+1} = 0, P~_L = 1, go DOWN with
//        P~_{l-1} = [ (2l+1) x P~_l - (l+1) P~_{l+1} ] / l,
//    then rescale with c = 1 / P~_0 so that P_0^(back) = 1.
static void backward_miller(double x, int lmax, int L, double *p)
{
    double *w = calloc((size_t)L + 2, sizeof(*w)); // w[L+1] = 0
    if (w == NULL) {
        fprintf(stderr, "Error: not enough memory\n");
        exit(1);
    }
    w[L] = 1.0;
    for (int l = L; l >= 1; l--)
        w[l - 1] = ((2.0 * l + 1.0) * x * w[l] - (l + 1.0) * w[l + 1]) / l;
    double c = 1.0 / w[0];
    for (int l = 0; l <= lmax; l++)
        p[l] = c * w[l];
    free(w);
}

// ****************************************************************************
// Helpers
// ****************************************************************************

static void usage(const char *prog)
{
    fprintf(stderr,
            "Usage: %s [l_max] [L] [x1 x2 ...]\n"
            "  l_max  largest degree l                                   (default %d)\n"
            "  L      starting degree of the backward recurrence, L > l_max (default %d)\n"
            "  x_i    points, with |x| < 1                (default 0.1 0.5 0.9 0.99)\n",
            prog, DEFAULT_L_MAX, DEFAULT_L);
}

static int parse_int(const char *s, int *out)
{
    char *end;
    errno = 0;
    long v = strtol(s, &end, 10);
    if (errno != 0 || end == s || *end != '\0' || v < 1 || v > 100000)
        return 0;
    *out = (int)v;
    return 1;
}

static int parse_double(const char *s, double *out)
{
    char *end;
    *out = strtod(s, &end);
    return end != s && *end == '\0';
}

// File-name tag for x: 0.99 -> "0p99", -0.5 -> "m0p5"
static void x_tag(double x, char *buf, size_t size)
{
    snprintf(buf, size, "%.6g", x);
    for (char *c = buf; *c; c++)
        if (*c == '.')
            *c = 'p';
        else if (*c == '-')
            *c = 'm';
}

// Relative error, or NAN when the reference is exactly zero.
static double rel_err(double approx, double ref)
{
    return ref != 0.0 ? (approx - ref) / ref : NAN;
}

// ****************************************************************************
// One x: compute, compare, write a CSV, print a summary line
// ****************************************************************************

static void study(double x, int lmax, int L)
{
    size_t n = (size_t)lmax + 1;
    double *fwd = malloc(n * sizeof(*fwd));
    double *bwd = malloc(n * sizeof(*bwd));
    double *gsl = malloc(n * sizeof(*gsl));
    quad *ref = malloc(n * sizeof(*ref));
    if (!fwd || !bwd || !gsl || !ref) {
        fprintf(stderr, "Error: not enough memory\n");
        exit(1);
    }

    quad xq = x; // the reference uses exactly the same (rounded) x as the double code
    forward_double(x, lmax, fwd);
    backward_miller(x, lmax, L, bwd);
    forward_quad(xq, lmax, ref);
    gsl_sf_legendre_Pl_array(lmax, x, gsl);

    char tag[32], path[128];
    x_tag(x, tag, sizeof(tag));
    snprintf(path, sizeof(path), OUTPUT_DIR "/legendre_P_x%s.csv", tag);
    FILE *f = fopen(path, "w");
    if (f == NULL) {
        perror(path);
        exit(1);
    }
    fprintf(f, "l,x,forward,backward_miller,reference_quad,gsl,abs_err_forward,rel_err_forward,"
               "abs_err_backward,rel_err_backward\n");

    double max_abs_f = 0, max_rel_f = 0, max_abs_b = 0, max_rel_b = 0;
    int gsl_equals_forward = 1;
    for (int l = 0; l <= lmax; l++) {
        double r = (double)ref[l];
        double ea_f = fabs(fwd[l] - r), er_f = rel_err(fwd[l], r);
        double ea_b = fabs(bwd[l] - r), er_b = rel_err(bwd[l], r);
        max_abs_f = fmax(max_abs_f, ea_f);
        max_abs_b = fmax(max_abs_b, ea_b);
        if (!isnan(er_f))
            max_rel_f = fmax(max_rel_f, fabs(er_f));
        if (!isnan(er_b))
            max_rel_b = fmax(max_rel_b, fabs(er_b));
        if (gsl[l] != fwd[l])
            gsl_equals_forward = 0;

        char refbuf[48];
        quadmath_snprintf(refbuf, sizeof(refbuf), "%.25Qe", ref[l]);
        fprintf(f, "%d,%.17g,%.17e,%.17e,%s,%.17e,%.6e,%.6e,%.6e,%.6e\n", l, x, fwd[l],
                bwd[l], refbuf, gsl[l], ea_f, er_f, ea_b, er_b);
    }
    fclose(f);

    char label[32];
    snprintf(label, sizeof(label), "P_l(%g)", x);
    printf("%-10s | %9.1e %9.1e | %9.1e %9.1e | %-3s | %s\n", label, max_abs_f, max_rel_f,
           max_abs_b, max_rel_b, gsl_equals_forward ? "yes" : "no", path);

    free(fwd);
    free(bwd);
    free(gsl);
    free(ref);
}

// ****************************************************************************

int main(int argc, char *argv[])
{
    int lmax = DEFAULT_L_MAX, L = DEFAULT_L;
    const double *xs = DEFAULT_X;
    size_t nx = sizeof(DEFAULT_X) / sizeof(DEFAULT_X[0]);
    double *user_x = NULL;

    if (argc >= 2 && !parse_int(argv[1], &lmax)) {
        fprintf(stderr, "Error: l_max must be a positive integer\n");
        usage(argv[0]);
        return 1;
    }
    if (argc >= 3 && !parse_int(argv[2], &L)) {
        fprintf(stderr, "Error: L must be a positive integer\n");
        usage(argv[0]);
        return 1;
    }
    if (lmax < 2 || L <= lmax) {
        fprintf(stderr, "Error: need l_max >= 2 and L > l_max\n");
        usage(argv[0]);
        return 1;
    }
    if (argc >= 4) {
        nx = (size_t)(argc - 3);
        user_x = malloc(nx * sizeof(*user_x));
        if (user_x == NULL)
            return 1;
        for (size_t i = 0; i < nx; i++)
            if (!parse_double(argv[3 + i], &user_x[i]) || fabs(user_x[i]) >= 1.0) {
                fprintf(stderr, "Error: x must be a number with |x| < 1, got '%s'\n",
                        argv[3 + i]);
                free(user_x);
                return 1;
            }
        xs = user_x;
    }

    gsl_set_error_handler_off(); // GSL errors are not fatal: we only compare its values

    if (mkdir(OUTPUT_DIR, 0775) != 0 && errno != EEXIST) {
        perror(OUTPUT_DIR);
        free(user_x);
        return 1;
    }

    printf("Legendre polynomials P_l(x): l = 0 ... %d, backward recurrence from L = %d\n", lmax,
           L);
    printf("Errors against the quad-precision reference (max over l)\n\n");
    printf("%-10s | %-19s | %-19s | %-3s |\n", "", "     forward", "  backward (Miller)", "GSL");
    printf("%-10s | %9s %9s | %9s %9s | %-3s | %s\n", "function", "abs err", "rel err",
           "abs err", "rel err", "=fw", "CSV file");
    printf("-----------+---------------------+---------------------+-----+---------\n");

    for (size_t i = 0; i < nx; i++)
        study(xs[i], lmax, L);

    printf("\n\"GSL =fw\": GSL gives exactly the same values as the forward recurrence.\n");
    free(user_x);
    return 0;
}
