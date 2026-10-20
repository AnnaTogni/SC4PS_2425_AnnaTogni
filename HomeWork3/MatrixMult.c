// HomeWork3
//
// Matrix multiplication  C = A B,   C_ij = sum_k A_ik B_kj   (i, j, k = 0 ... N-1)
//
//   A = a * 1_N ,  B = b * 1_N
//   1_N = N x N matrix of ones (default)      -> every C_ij = N a b
//       or N x N identity matrix ("identity") -> C = a b I
//
// The product is computed with all 6 loop orderings (ijk, ikj, jik, jki, kij, kji),
// each one is timed (benchmark), every result is checked, and C is saved to a file.
//
// ******************************************************
// HOW TO COMPILE:
//   make                     (or: gcc -Wall -Wextra -O3 MatrixMult.c -o matrixmult -lm)
//
// HOW TO RUN:
//   ./matrixmult a b N fileout [ones|identity]
//   ./matrixmult 2 3 1000 C.dat
//   ./matrixmult 2 3 1000 C.dat identity
// ******************************************************

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Matrices are stored "row-major" in a 1D array: element (i, j) is M[i*N + j].
// Consecutive j are next to each other in memory, consecutive i are N doubles apart.
#define IDX(i, j) ((i) * N + (j))

// Largest matrix written to the terminal (bigger ones only go to the file).
#define N_PRINT 6

static void usage(const char *prog)
{
    fprintf(stderr,
            "Usage: %s a b N fileout [ones|identity]\n"
            "  a, b      value of the elements of A and B\n"
            "  N         matrix dimension (N x N)\n"
            "  fileout   output file where C is saved\n"
            "  ones      (default) A, B have every element equal to a, b -> C_ij = N a b\n"
            "  identity  A = a I, B = b I                                -> C = a b I\n"
            "Example: %s 2 3 1000 C.dat\n",
            prog, prog);
}

static int parse_double(const char *s, double *out)
{
    char *end;
    *out = strtod(s, &end);
    return end != s && *end == '\0';
}

static double seconds_now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

// ****************************************************************************
// The 6 loop orderings. They do exactly the same arithmetic, only the order in
// which the memory is visited changes. C must be zero before the call.
// ****************************************************************************

static void matmul_ijk(size_t N, const double *A, const double *B, double *C)
{
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++)
            for (size_t k = 0; k < N; k++)
                C[IDX(i, j)] += A[IDX(i, k)] * B[IDX(k, j)];
}

static void matmul_ikj(size_t N, const double *A, const double *B, double *C)
{
    for (size_t i = 0; i < N; i++)
        for (size_t k = 0; k < N; k++)
            for (size_t j = 0; j < N; j++)
                C[IDX(i, j)] += A[IDX(i, k)] * B[IDX(k, j)];
}

static void matmul_jik(size_t N, const double *A, const double *B, double *C)
{
    for (size_t j = 0; j < N; j++)
        for (size_t i = 0; i < N; i++)
            for (size_t k = 0; k < N; k++)
                C[IDX(i, j)] += A[IDX(i, k)] * B[IDX(k, j)];
}

static void matmul_jki(size_t N, const double *A, const double *B, double *C)
{
    for (size_t j = 0; j < N; j++)
        for (size_t k = 0; k < N; k++)
            for (size_t i = 0; i < N; i++)
                C[IDX(i, j)] += A[IDX(i, k)] * B[IDX(k, j)];
}

static void matmul_kij(size_t N, const double *A, const double *B, double *C)
{
    for (size_t k = 0; k < N; k++)
        for (size_t i = 0; i < N; i++)
            for (size_t j = 0; j < N; j++)
                C[IDX(i, j)] += A[IDX(i, k)] * B[IDX(k, j)];
}

static void matmul_kji(size_t N, const double *A, const double *B, double *C)
{
    for (size_t k = 0; k < N; k++)
        for (size_t j = 0; j < N; j++)
            for (size_t i = 0; i < N; i++)
                C[IDX(i, j)] += A[IDX(i, k)] * B[IDX(k, j)];
}

typedef void (*matmul_fn)(size_t, const double *, const double *, double *);

static const struct {
    const char *name;
    matmul_fn   fn;
    const char *inner; // what the innermost loop walks through
} ORDERINGS[] = {
    {"ijk", matmul_ijk, "A row, B column"},
    {"ikj", matmul_ikj, "B row, C row    (all contiguous)"},
    {"jik", matmul_jik, "A row, B column"},
    {"jki", matmul_jki, "A column, C column (all strided)"},
    {"kij", matmul_kij, "B row, C row    (all contiguous)"},
    {"kji", matmul_kji, "A column, C column (all strided)"},
};
#define N_ORDERINGS (sizeof(ORDERINGS) / sizeof(ORDERINGS[0]))

// ****************************************************************************
// Checks
// ****************************************************************************

// Value that element (i, j) of C must have.
static double expected_value(size_t i, size_t j, size_t N, double a, double b, int identity)
{
    if (identity)
        return (i == j) ? a * b : 0.0;
    return (double)N * a * b;
}

// Full check: compare EVERY element with the expected value, within a tolerance.
// Each C_ij is a sum of N terms, and every addition can add a rounding error of
// about DBL_EPSILON, so the error can grow up to ~ N * DBL_EPSILON * |C_ij|.
// Returns the number of wrong elements.
static size_t check_all(const double *C, size_t N, double a, double b, int identity,
                        size_t *n_exact)
{
    double tol = (double)N * DBL_EPSILON * fabs((double)N * a * b);
    size_t n_wrong = 0;
    *n_exact = 0;
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++) {
            double e = expected_value(i, j, N, a, b, identity);
            double c = C[IDX(i, j)];
            if (c == e)
                (*n_exact)++;
            if (fabs(c - e) > tol)
                n_wrong++;
        }
    return n_wrong;
}

// Quick check: a single comparison instead of N^2 of them.
// The sum of all the elements of C must be
//   ones:     N^2 elements * N a b = N^3 a b
//   identity: N diagonal elements * a b = N a b
// It costs one pass over C but only ONE comparison, and needs no knowledge of
// the single elements. (It can miss errors that cancel out, e.g. +1 and -1.)
static int quick_check(const double *C, size_t N, double a, double b, int identity,
                       double *sum_out, double *expected_out)
{
    double sum = 0.0;
    for (size_t n = 0; n < N * N; n++)
        sum += C[n];
    double expected = identity ? (double)N * a * b : (double)N * N * N * a * b;
    *sum_out = sum;
    *expected_out = expected;
    return fabs(sum - expected) <= (double)N * N * DBL_EPSILON * fabs(expected);
}

// ****************************************************************************

static int save_matrix(const char *filename, const double *C, size_t N)
{
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        perror(filename);
        return 0;
    }
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++)
            fprintf(f, "%.17g%c", C[IDX(i, j)], j + 1 < N ? ' ' : '\n');
    }
    if (fclose(f) != 0) {
        perror(filename);
        return 0;
    }
    return 1;
}

static void print_matrix(const char *name, const double *M, size_t N)
{
    printf("%s =\n", name);
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++)
            printf(" %8.3g", M[IDX(i, j)]);
        printf("\n");
    }
}

int main(int argc, char *argv[])
{
    // ---- 1. Input from the command line: a b N fileout [ones|identity] -----
    if (argc != 5 && argc != 6) {
        usage(argv[0]);
        return 1;
    }

    double a, b, N_in;
    if (!parse_double(argv[1], &a) || !parse_double(argv[2], &b)) {
        fprintf(stderr, "Error: a and b must be numbers\n");
        usage(argv[0]);
        return 1;
    }
    if (!parse_double(argv[3], &N_in) || N_in < 1 || N_in != floor(N_in)) {
        fprintf(stderr, "Error: N must be a positive integer, got '%s'\n", argv[3]);
        return 1;
    }
    const char *fileout = argv[4];

    int identity = 0;
    if (argc == 6) {
        if (strcmp(argv[5], "identity") == 0)
            identity = 1;
        else if (strcmp(argv[5], "ones") != 0) {
            fprintf(stderr, "Error: last argument must be 'ones' or 'identity'\n");
            return 1;
        }
    }
    size_t N = (size_t)N_in;

    printf("C = A B  with  A = %g * %s,  B = %g * %s,  N = %zu\n", a,
           identity ? "I" : "ones", b, identity ? "I" : "ones", N);
    printf("Memory: 3 matrices x %zu x %zu x 8 bytes = %.1f MB\n", N, N,
           3.0 * N * N * sizeof(double) / 1e6);
    printf("Operations: 2 N^3 = %.3g floating-point operations per product\n\n",
           2.0 * N * N * N);
    fflush(stdout);

    // ---- 2. Allocate and fill A and B ---------------------------------------
    double *A = malloc(N * N * sizeof(double));
    double *B = malloc(N * N * sizeof(double));
    double *C = malloc(N * N * sizeof(double));
    if (A == NULL || B == NULL || C == NULL) {
        fprintf(stderr, "Error: not enough memory for N = %zu\n", N);
        free(A);
        free(B);
        free(C);
        return 1;
    }
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++) {
            A[IDX(i, j)] = (identity && i != j) ? 0.0 : a;
            B[IDX(i, j)] = (identity && i != j) ? 0.0 : b;
        }

    // ---- 3. Benchmark the 6 loop orderings, checking each result -----------
    printf("Ordering  time [s]   GFLOP/s   all C_ij correct?   inner loop walks along\n");
    printf("--------  --------   -------   -----------------   ----------------------\n");

    double times[N_ORDERINGS];
    int all_ok = 1;
    size_t best = 0;
    size_t n_exact = 0;
    for (size_t o = 0; o < N_ORDERINGS; o++) {
        memset(C, 0, N * N * sizeof(double));

        double t0 = seconds_now();
        ORDERINGS[o].fn(N, A, B, C);
        times[o] = seconds_now() - t0;

        size_t n_wrong = check_all(C, N, a, b, identity, &n_exact);
        all_ok = all_ok && (n_wrong == 0);
        if (times[o] < times[best])
            best = o;

        printf("  %s    %8.4f   %7.2f   %-17s   %s\n", ORDERINGS[o].name, times[o],
               2.0 * N * N * N / times[o] / 1e9, n_wrong == 0 ? "YES" : "NO",
               ORDERINGS[o].inner);
        fflush(stdout);
    }

    size_t worst = 0;
    for (size_t o = 1; o < N_ORDERINGS; o++)
        if (times[o] > times[worst])
            worst = o;
    printf("\nFastest: %s (%.4f s)   slowest: %s (%.4f s)   ratio: %.1fx\n",
           ORDERINGS[best].name, times[best], ORDERINGS[worst].name, times[worst],
           times[worst] / times[best]);

    // ---- 4. Checks on the final C ------------------------------------------
    double e = expected_value(0, 0, N, a, b, identity);
    printf("\nFull check  (all %zu elements): %s\n", N * N, all_ok ? "PASSED" : "FAILED");
    printf("  expected C_00 = %.17g, computed C_00 = %.17g\n", e, C[0]);
    printf("  elements exactly equal (==) to the expected value: %zu of %zu\n",
           n_exact, N * N);

    double sum, sum_expected;
    int quick_ok = quick_check(C, N, a, b, identity, &sum, &sum_expected);
    printf("Quick check (1 comparison, sum of all C_ij = %s): %s\n",
           identity ? "N a b" : "N^3 a b", quick_ok ? "PASSED" : "FAILED");
    printf("  sum = %.17g, expected = %.17g\n", sum, sum_expected);

    if (N <= N_PRINT) {
        printf("\n");
        print_matrix("A", A, N);
        print_matrix("B", B, N);
        print_matrix("C", C, N);
    }

    // ---- 5. Save C ------------------------------------------------------------
    int saved = save_matrix(fileout, C, N);
    if (saved)
        printf("\nMatrix C saved to '%s'\n", fileout);

    free(A);
    free(B);
    free(C);

    if (!saved)
        return 1;
    return (all_ok && quick_ok) ? 0 : 2;
}
