// legendre_stability.c
// This program computes Legendre polynomials P_l(x) for a range of l values and specific x values
// using both a forward recursion and Miller's backward method.
// It compares the results to reference values computed with the GSL special functions library.
// The results are written to CSV files for further analysis.

#include <float.h>
#include <errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <math.h>
#include <quadmath.h>
#include <stdio.h>
#include <stdlib.h>

#define DEFAULT_L_MAX 50
#define DEFAULT_MILLER_L 100
#define DEFAULT_NUM_X_VALUES 4

static const double DEFAULT_X_VALUES[DEFAULT_NUM_X_VALUES] = {0.1, 0.5, 0.9, 0.99};

// Utility functions for argument parsing, error calculation, and CSV output formatting.

// Returns the absolute value of a double.
static double abs_double(double x) {
    return x < 0.0 ? -x : x;
}

// Parses a string as a non-negative integer with error checking.
static int parse_int_arg(const char *text, const char *name) {
    char *endptr = NULL;
    long value = 0;

    errno = 0;
    value = strtol(text, &endptr, 10);
    if (errno != 0 || endptr == text || *endptr != '\0' || value < 0 || value > 1000000L) {
        fprintf(stderr, "Invalid %s: %s\n", name, text);
        exit(EXIT_FAILURE);
    }

    return (int)value;
}

// Parses a string as a double with error checking.
static double parse_double_arg(const char *text, const char *name) {
    char *endptr = NULL;
    double value = 0.0;

    errno = 0;
    value = strtod(text, &endptr);
    if (errno != 0 || endptr == text || *endptr != '\0') {
        fprintf(stderr, "Invalid %s: %s\n", name, text);
        exit(EXIT_FAILURE);
    }

    return value;
}

// Converts a double x into a string tag suitable for filenames by replacing '.' with 'p', '-' with 'm', and '+' with 'p'.
static void make_x_tag(double x, char *buffer, size_t buffer_size) {
    int written = snprintf(buffer, buffer_size, "%.6g", x);
    if (written < 0 || (size_t)written >= buffer_size) {
        fprintf(stderr, "Failed to format x value\n");
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; buffer[i] != '\0'; ++i) {
        if (buffer[i] == '.') {
            buffer[i] = 'p';
        } else if (buffer[i] == '-') {
            buffer[i] = 'm';
        } else if (buffer[i] == '+') {
            buffer[i] = 'p';
        }
    }
}

// Computes the signed error between an approximation and a reference value.
static double error_value(double approx, double reference) {
    return approx - reference;
}

// Computes the signed relative error between an approximation and a reference value,
// returning NAN when the reference scale is too small for a meaningful relative error.
static double rel_error(double approx, double reference) {
    if (abs_double(reference) < DBL_MIN) {
        return NAN;
    }
    return error_value(approx, reference) / reference;
}

// Computes Legendre polynomials P_l(x) for l = 0 to lmax using the forward recursion relation.
static void legendre_forward_double(double x, int lmax, double *values) {
    values[0] = 1.0;
    if (lmax == 0) {
        return;
    }
    values[1] = x;
    for (int ell = 1; ell < lmax; ++ell) {
        values[ell + 1] = ((2.0 * ell + 1.0) * x * values[ell] - ell * values[ell - 1]) / (ell + 1.0);
    }
}

// Computes reference values P_l(x) for l = 0 to lmax using the GSL special functions library.
static void legendre_reference_gsl(double x, int lmax, double *values) {
    int status = gsl_sf_legendre_Pl_array(lmax, x, values);
    if (status != 0) {
        fprintf(stderr, "GSL failed to compute Legendre reference values for x = %.17g\n", x);
        exit(EXIT_FAILURE);
    }
}

// Computes reference values P_l(x) for l = 0 to lmax in quad precision.
static void legendre_reference_quad(double x, int lmax, __float128 *values) {
    __float128 xq = (__float128)x;

    values[0] = (__float128)1.0;
    if (lmax == 0) {
        return;
    }

    values[1] = xq;
    for (int ell = 1; ell < lmax; ++ell) {
        values[ell + 1] =
            (((__float128)(2 * ell + 1)) * xq * values[ell] - ((__float128)ell) * values[ell - 1]) /
            ((__float128)(ell + 1));
    }
}

// Computes Legendre polynomials P_l(x) for l = 0 to lmax using Miller's backward method starting from start_l > lmax.
static void legendre_backward_miller(double x, int lmax, int start_l, double *values) {
    double *work = calloc((size_t)start_l + 2U, sizeof(*work));
    if (work == NULL) {
        fprintf(stderr, "Allocation failed in legendre_backward_miller\n");
        exit(EXIT_FAILURE);
    }

    work[start_l + 1] = 0.0;
    work[start_l] = 1.0;

    for (int ell = start_l; ell >= 1; --ell) {
        work[ell - 1] = ((2.0 * ell + 1.0) * x * work[ell] - (ell + 1.0) * work[ell + 1]) / ell;
    }

    {
        double scale = 1.0 / work[0];
        for (int ell = 0; ell <= lmax; ++ell) {
            values[ell] = scale * work[ell];
        }
    }

    free(work);
}

// For a given x, computes Legendre polynomials using both methods, compares them to the quad/GSL references,
// and writes the results to a CSV file.
static void write_csv_for_x(double x, int lmax, int start_l) {

    // Prepare filename and allocate arrays for results
    char filename[64];
    char x_tag[32];
    FILE *output = NULL;
    double *forward = calloc((size_t)lmax + 1U, sizeof(*forward));
    double *backward = calloc((size_t)lmax + 1U, sizeof(*backward));
    double *reference_gsl = calloc((size_t)lmax + 1U, sizeof(*reference_gsl));
    __float128 *reference_quad = calloc((size_t)lmax + 1U, sizeof(*reference_quad));

    if (forward == NULL || backward == NULL || reference_gsl == NULL || reference_quad == NULL) {
        fprintf(stderr, "Allocation failed in write_csv_for_x\n");
        free(forward);
        free(backward);
        free(reference_gsl);
        free(reference_quad);
        exit(EXIT_FAILURE);
    }

    // Create a filename based on the x value and open the CSV file for writing
    make_x_tag(x, x_tag, sizeof(x_tag));
    snprintf(filename, sizeof(filename), "legendre_x_%s.csv", x_tag);
    output = fopen(filename, "w");
    if (output == NULL) {
        fprintf(stderr, "Cannot open %s for writing\n", filename);
        free(forward);
        free(backward);
        free(reference_gsl);
        free(reference_quad);
        exit(EXIT_FAILURE);
    }

    // Compute Legendre polynomials using both methods and the reference values
    legendre_forward_double(x, lmax, forward);
    legendre_backward_miller(x, lmax, start_l, backward);
    legendre_reference_gsl(x, lmax, reference_gsl);
    legendre_reference_quad(x, lmax, reference_quad);

    // Write CSV header
    fprintf(output,
            "l,x,forward_double,backward_miller,reference_quad,reference_gsl,err_forward,rel_err_forward,"
            "err_backward,rel_err_backward\n");

    // Write results for each l to the CSV file, including signed and relative errors
    for (int ell = 0; ell <= lmax; ++ell) {

        // Calculate errors compared to the quad-precision reference
        // double ref = (double)reference_quad[ell];
        double ref = (double)reference_gsl[ell];
        double ref_quad = (double)reference_quad[ell];
        double err_f = error_value(forward[ell], ref);
        double rel_err_f = rel_error(forward[ell], ref);
        double err_b = error_value(backward[ell], ref);
        double rel_err_b = rel_error(backward[ell], ref);

        fprintf(output,
                "%d,%.2f,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
                ell,
                x,
                forward[ell],
                backward[ell],
                ref_quad,
                reference_gsl[ell],
                err_f,
                rel_err_f,
                err_b,
                rel_err_b);
    }

    fclose(output);
    free(forward);
    free(backward);
    free(reference_gsl);
    free(reference_quad);
}

// For a given x, computes Legendre polynomials using both methods, compares them to the quad reference,
// and prints a summary line with the maximum relative errors for both methods.
static void print_summary_line(double x, int lmax, int start_l) {

    // Prepare arrays for results and variables to track maximum relative errors
    char x_tag[32];
    double *forward = calloc((size_t)lmax + 1U, sizeof(*forward));
    double *backward = calloc((size_t)lmax + 1U, sizeof(*backward));
    __float128 *reference_quad = calloc((size_t)lmax + 1U, sizeof(*reference_quad));
    double max_rel_forward = 0.0;
    double max_rel_backward = 0.0;

    if (forward == NULL || backward == NULL || reference_quad == NULL) {
        fprintf(stderr, "Allocation failed in print_summary_line\n");
        free(forward);
        free(backward);
        free(reference_quad);
        exit(EXIT_FAILURE);
    }

    // Compute Legendre polynomials using both methods and the quad reference values
    legendre_forward_double(x, lmax, forward);
    legendre_backward_miller(x, lmax, start_l, backward);
    legendre_reference_quad(x, lmax, reference_quad);
    make_x_tag(x, x_tag, sizeof(x_tag));

    // Calculate the maximum absolute relative error for both methods across all l values
    for (int ell = 0; ell <= lmax; ++ell) {
        double ref = (double)reference_quad[ell];
        double rel_f = rel_error(forward[ell], ref);
        double rel_b = rel_error(backward[ell], ref);
        if (!isnan(rel_f) && abs_double(rel_f) > max_rel_forward) {
            max_rel_forward = abs_double(rel_f);
        }
        if (!isnan(rel_b) && abs_double(rel_b) > max_rel_backward) {
            max_rel_backward = abs_double(rel_b);
        }
    }

    // Print a summary line with the maximum relative errors for both methods
    printf("x = %.6g | max rel. err forward = %.3e | max rel. err Miller = %.3e | CSV: legendre_x_%s.csv\n",
           x,
           max_rel_forward,
           max_rel_backward,
           x_tag);

    free(forward);
    free(backward);
    free(reference_quad);
}

// The main function parses command-line arguments for l_max, Miller start L, and x values, 
// then computes and outputs the stability study results for each x value.
int main(int argc, char **argv) {

    // Default parameters for l_max, Miller start L, and x values
    int lmax = DEFAULT_L_MAX;
    int miller_l = DEFAULT_MILLER_L;
    const double *x_values = DEFAULT_X_VALUES;
    int num_x_values = DEFAULT_NUM_X_VALUES;
    double *user_x_values = NULL;

    // Parse command-line arguments to override defaults for l_max, Miller start L, and x values
    if (argc >= 2) {
        lmax = parse_int_arg(argv[1], "l_max");
    }
    if (argc >= 3) {
        miller_l = parse_int_arg(argv[2], "Miller start L");
    }
    if (argc >= 4) {
        num_x_values = argc - 3;
        user_x_values = calloc((size_t)num_x_values, sizeof(*user_x_values));
        if (user_x_values == NULL) {
            fprintf(stderr, "Allocation failed while reading x values\n");
            return EXIT_FAILURE;
        }
        for (int i = 0; i < num_x_values; ++i) {
            user_x_values[i] = parse_double_arg(argv[i + 3], "x");
        }
        x_values = user_x_values;
    }
    if (miller_l <= lmax) {
        fprintf(stderr, "Miller start L must be larger than l_max\n");
        free(user_x_values);
        return EXIT_FAILURE;
    }

    printf("Legendre polynomial stability study\n");
    printf("Parameters: l_max = %d, Miller start L = %d\n\n", lmax, miller_l);

    // For each x value, compute the Legendre polynomials using both methods, compare to GSL reference, 
    // and write results to CSV files.
    for (int i = 0; i < num_x_values; ++i) {
        write_csv_for_x(x_values[i], lmax, miller_l);
        print_summary_line(x_values[i], lmax, miller_l);
    }

    printf("\nEach CSV contains P_l(x), the quad/GSL references, and signed/relative errors.\n");
    printf("You can plot err_forward, err_backward, rel_err_forward, and rel_err_backward versus l for each x.\n");
    printf("Usage: ./legendre_stability [l_max] [miller_L] [x1 x2 ...]\n");

    free(user_x_values);
    return 0;
}
