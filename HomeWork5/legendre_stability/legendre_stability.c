// HomeWork5 -- Stability of the three-term recurrence for Legendre functions
//
// The recurrence (Bonnet)
//     (l+1) f_{l+1}(x) = (2l+1) x f_l(x) - l f_{l-1}(x)
// has two independent solutions: P_l(x) (first kind) and Q_l(x) (second kind).
// Which one is computed depends only on the two starting values.
//
// For every x the program compares, for l = 0 ... l_max:
//   - forward recurrence (l increasing) in double precision
//   - Miller's backward recurrence (l decreasing, from L > l_max) in double,
//     normalized with ONE condition (f_0), as in the HomeWork05 text
//   - for |x| < 1 only: backward recurrence normalized with TWO conditions
//     (P_0 = 1 and P_1 = x), as in the professor's corrected example (lecture 8)
// against a reference computed in quad precision (__float128, 113-bit mantissa).
//
//   |x| < 1 : the function studied is P_l(x)   (the homework's case)
//             reference = forward recurrence in quad precision
//   x  > 1 : the function studied is Q_l(x)   (where the recurrence really is unstable)
//             reference = Miller backward recurrence in quad precision from l_max + 1000
//
// GSL (gsl_sf_legendre_Pl_array / gsl_sf_legendre_Ql) is also printed as a cross-check.
//
// ******************************************************
// HOW TO COMPILE:
//   make
// HOW TO RUN:
//   ./legendre_stability [l_max] [L] [x1 x2 ...]
//   ./legendre_stability                  (l_max = 50, L = 80, x = 0.1 0.5 0.9 0.99 1.1 1.5 3)
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
#define REF_EXTRA     1000 // the quad Miller reference for Q_l starts at l_max + REF_EXTRA

static const double DEFAULT_X[] = {0.1, 0.5, 0.9, 0.99, 1.1, 1.5, 3.0};

typedef __float128 quad;

// ****************************************************************************
// Recurrences
// ****************************************************************************

// Forward recurrence from the starting values f[0], f[1] (already set), up to f[lmax].
static void forward_double(double x, int lmax, double *f)
{
    for (int l = 1; l < lmax; l++)
        f[l + 1] = ((2.0 * l + 1.0) * x * f[l] - l * f[l - 1]) / (l + 1.0);
}

static void forward_quad(quad x, int lmax, quad *f)
{
    for (int l = 1; l < lmax; l++)
        f[l + 1] = ((quad)(2 * l + 1) * x * f[l] - (quad)l * f[l - 1]) / (quad)(l + 1);
}

// Miller: start from f[L+1] = 0, f[L] = 1, go DOWN with
//     f_{l-1} = [ (2l+1) x f_l - (l+1) f_{l+1} ] / l,
// then rescale everything so that f[0] = f0 (the known value of the function at l = 0).
static void miller_double(double x, int lmax, int L, double f0, double *f)
{
    double *w = calloc((size_t)L + 2, sizeof(*w));
    if (w == NULL) {
        fprintf(stderr, "Error: not enough memory\n");
        exit(1);
    }
    w[L] = 1.0;
    for (int l = L; l >= 1; l--)
        w[l - 1] = ((2.0 * l + 1.0) * x * w[l] - (l + 1.0) * w[l + 1]) / l;
    double c = f0 / w[0];
    for (int l = 0; l <= lmax; l++)
        f[l] = c * w[l];
    free(w);
}

// Backward recurrence with TWO normalization conditions f[0] = f0, f[1] = f1.
// Every solution of the recurrence is a P_l + b Q_l. Two backward runs from independent
// starts (f_L = 1, f_{L+1} = 0 and f_L = 0, f_{L+1} = 1) span all solutions; the
// combination c_a A + c_b B with the right f_0 and f_1 is the wanted one.
// Works when the two solutions have similar size (|x| < 1); for x > 1 the two runs become
// proportional (both ~ Q_l) and the 2x2 system is singular.
static void backward_two_conditions_double(double x, int lmax, int L, double f0, double f1,
                                           double *f)
{
    double *a = calloc((size_t)L + 2, sizeof(*a));
    double *b = calloc((size_t)L + 2, sizeof(*b));
    if (a == NULL || b == NULL) {
        fprintf(stderr, "Error: not enough memory\n");
        exit(1);
    }
    a[L] = 1.0; // start A: f_L = 1, f_{L+1} = 0
    b[L + 1] = 1.0; // start B: f_L = 0, f_{L+1} = 1
    for (int l = L; l >= 1; l--) {
        a[l - 1] = ((2.0 * l + 1.0) * x * a[l] - (l + 1.0) * a[l + 1]) / l;
        b[l - 1] = ((2.0 * l + 1.0) * x * b[l] - (l + 1.0) * b[l + 1]) / l;
    }
    double det = a[0] * b[1] - a[1] * b[0];
    double ca = (f0 * b[1] - f1 * b[0]) / det;
    double cb = (f1 * a[0] - f0 * a[1]) / det;
    for (int l = 0; l <= lmax; l++)
        f[l] = ca * a[l] + cb * b[l];
    free(a);
    free(b);
}

static void miller_quad(quad x, int lmax, int L, quad f0, quad *f)
{
    quad *w = calloc((size_t)L + 2, sizeof(*w));
    if (w == NULL) {
        fprintf(stderr, "Error: not enough memory\n");
        exit(1);
    }
    w[L] = 1;
    for (int l = L; l >= 1; l--)
        w[l - 1] = ((quad)(2 * l + 1) * x * w[l] - (quad)(l + 1) * w[l + 1]) / (quad)l;
    quad c = f0 / w[0];
    for (int l = 0; l <= lmax; l++)
        f[l] = c * w[l];
    free(w);
}

// ****************************************************************************
// Helpers
// ****************************************************************************

static void usage(const char *prog)
{
    fprintf(stderr,
            "Usage: %s [l_max] [L] [x1 x2 ...]\n"
            "  l_max  largest degree l              (default %d)\n"
            "  L      starting degree of Miller's backward recurrence, L > l_max (default %d)\n"
            "  x_i    points: |x| < 1 studies P_l(x), x > 1 studies Q_l(x)\n"
            "         (default 0.1 0.5 0.9 0.99 1.1 1.5 3)\n",
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
// One study: compute, compare, write a CSV, print a summary line
// ****************************************************************************

static int study(double x, int lmax, int L)
{
    int is_P = fabs(x) < 1.0;
    const char *name = is_P ? "P" : "Q";
    size_t n = (size_t)lmax + 1;

    double *fwd = malloc(n * sizeof(*fwd));
    double *bwd = malloc(n * sizeof(*bwd));
    double *bw2 = malloc(n * sizeof(*bw2));
    double *gsl = malloc(n * sizeof(*gsl));
    double *cond = malloc(n * sizeof(*cond));
    quad *ref = malloc(n * sizeof(*ref));
    if (!fwd || !bwd || !bw2 || !gsl || !cond || !ref) {
        fprintf(stderr, "Error: not enough memory\n");
        exit(1);
    }

    quad xq = x; // the reference uses exactly the same (rounded) x as the double code
    if (is_P) {
        // P_0 = 1, P_1 = x
        fwd[0] = 1.0;
        fwd[1] = x;
        forward_double(x, lmax, fwd);
        miller_double(x, lmax, L, 1.0, bwd);
        backward_two_conditions_double(x, lmax, L, 1.0, x, bw2);
        ref[0] = 1;
        ref[1] = xq;
        forward_quad(xq, lmax, ref);
        gsl_sf_legendre_Pl_array(lmax, x, gsl);
    } else {
        // Q_0 = 1/2 ln((x+1)/(x-1)), Q_1 = x Q_0 - 1
        double q0 = 0.5 * log((x + 1.0) / (x - 1.0));
        quad q0q = (quad)0.5 * logq((xq + 1) / (xq - 1));
        fwd[0] = q0;
        fwd[1] = x * q0 - 1.0;
        forward_double(x, lmax, fwd);
        miller_double(x, lmax, L, q0, bwd);
        for (int l = 0; l <= lmax; l++)
            bw2[l] = NAN; // two conditions: singular for x > 1 (see the function)
        miller_quad(xq, lmax, lmax + REF_EXTRA, q0q, ref);
        for (int l = 0; l <= lmax; l++)
            gsl[l] = gsl_sf_legendre_Ql(l, x);
    }

    // Relative condition number of P_l at x:  kappa = |x P_l'(x) / P_l(x)|,
    // with P_l'(x) = l (x P_l - P_{l-1}) / (x^2 - 1). It says how much a relative
    // error in x (e.g. from rounding x to a double) is amplified in P_l(x).
    for (int l = 0; l <= lmax; l++) {
        if (!is_P || l == 0) {
            cond[l] = is_P ? 0.0 : NAN;
            continue;
        }
        quad dP = (quad)l * (xq * ref[l] - ref[l - 1]) / (xq * xq - 1);
        cond[l] = ref[l] != 0 ? (double)fabsq(xq * dP / ref[l]) : INFINITY;
    }

    char tag[32], path[128];
    x_tag(x, tag, sizeof(tag));
    snprintf(path, sizeof(path), OUTPUT_DIR "/legendre_%s_x%s.csv", name, tag);
    FILE *f = fopen(path, "w");
    if (f == NULL) {
        perror(path);
        exit(1);
    }
    fprintf(f, "l,x,forward,backward_miller,reference_quad,gsl,abs_err_forward,rel_err_forward,"
               "abs_err_backward,rel_err_backward,condition_number,backward_two_conditions,"
               "abs_err_backward2,rel_err_backward2\n");

    double max_abs_f = 0, max_rel_f = 0, max_abs_b = 0, max_rel_b = 0, max_abs_2 = 0;
    int gsl_equals_forward = 1;
    for (int l = 0; l <= lmax; l++) {
        double r = (double)ref[l];
        double ea_f = fabs(fwd[l] - r), er_f = rel_err(fwd[l], r);
        double ea_b = fabs(bwd[l] - r), er_b = rel_err(bwd[l], r);
        double ea_2 = fabs(bw2[l] - r), er_2 = rel_err(bw2[l], r);
        if (is_P)
            max_abs_2 = fmax(max_abs_2, ea_2);
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
        fprintf(f, "%d,%.17g,%.17e,%.17e,%s,%.17e,%.6e,%.6e,%.6e,%.6e,%.6e,%.17e,%.6e,%.6e\n", l,
                x, fwd[l], bwd[l], refbuf, gsl[l], ea_f, er_f, ea_b, er_b, cond[l], bw2[l],
                ea_2, er_2);
    }
    fclose(f);

    char label[32];
    snprintf(label, sizeof(label), "%s_l(%g)", name, x);
    char two[16];
    if (is_P)
        snprintf(two, sizeof(two), "%9.1e", max_abs_2);
    else
        snprintf(two, sizeof(two), "%9s", "-");
    printf("%-10s | %9.1e %9.1e | %9.1e %9.1e | %s | %-3s | %s\n", label, max_abs_f, max_rel_f,
           max_abs_b, max_rel_b, two, gsl_equals_forward ? "yes" : "no", path);

    free(fwd);
    free(bwd);
    free(bw2);
    free(gsl);
    free(cond);
    free(ref);
    return 0;
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
            if (!parse_double(argv[3 + i], &user_x[i]) || user_x[i] <= -1.0 ||
                user_x[i] == 1.0) {
                fprintf(stderr, "Error: x must be a number with x > -1 and x != 1, got '%s'\n",
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

    printf("Legendre recurrence stability: l = 0 ... %d, Miller start L = %d\n", lmax, L);
    printf("Errors against the quad-precision reference (max over l)\n\n");
    printf("%-10s | %-19s | %-19s | %-9s | %-3s |\n", "", "     forward", " backward, f_0 only",
           " bw P0+P1", "GSL");
    printf("%-10s | %9s %9s | %9s %9s | %9s | %-3s | %s\n", "function", "abs err", "rel err",
           "abs err", "rel err", "abs err", "=fw", "CSV file");
    printf("-----------+---------------------+---------------------+-----------+-----+---------\n");

    for (size_t i = 0; i < nx; i++)
        study(xs[i], lmax, L);

    printf("\n|x| < 1: P_l(x).  x > 1: Q_l(x).  \"GSL =fw\": GSL gives exactly the forward values.\n");
    printf("\"backward, f_0 only\" = Miller as in the HomeWork05 text; \"bw P0+P1\" = backward\n"
           "recurrence normalized with P_0 = 1 and P_1 = x, as in the professor's correction.\n");
    free(user_x);
    return 0;
}
