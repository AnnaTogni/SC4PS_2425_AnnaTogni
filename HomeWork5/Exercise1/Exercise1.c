// 
// Exercise 1: Summation of floating-point numbers
// This program demonstrates the effect of floating-point addition order on the result.
// It compares three methods of summation: 
// a simple for loop, the Kahan summation algorithm, and the GSL library function (if available).
// 
// The expected result is 0.5, but due to the nature of floating-point arithmetic, 
// the results may differ based on the method used.
// 
// Compile with: gcc -o Exercise1 Exercise1.c -lm
// To enable GSL summation, compile with: gcc -o Exercise1 Exercise1.c -lm -lgsl -lgslcblas -DUSE_GSL
// 


#include <stdio.h>
#include <stddef.h>

#ifdef USE_GSL
#include <gsl/gsl_vector.h>
#endif
// This means:
// - only include the GSL header if USE_GSL is defined at compile time
// - if you compile normally, this part is ignored
// - if you compile with something like -DUSE_GSL, then the GSL code becomes active
// So this file works both with and without GSL installed.

#define SUM_FOR(result, arr, n)                  \
    do {                                         \
        (result) = 0.0;                          \
        for (size_t i = 0; i < (n); ++i) {       \
            (result) += (arr)[i];                \
        }                                        \
    } while (0)

#define SUM_KAHAN(result, arr, n)                \
    do {                                         \
        double sum_ = 0.0;                       \
        double c_ = 0.0;                         \
        for (size_t i = 0; i < (n); ++i) {       \
            double y_ = (arr)[i] - c_;           \
            double t_ = sum_ + y_;               \
            c_ = (t_ - sum_) - y_;               \
            sum_ = t_;                           \
        }                                        \
        (result) = sum_;                         \
    } while (0)

// Both SUM_FOR and SUM_KAHAN are C preprocessor macros, not functions.
// That means before compilation, the preprocessor literally replaces the macro call with its body.
// Each line ends with \ so the macro continues onto the next line. 
// Without those backslashes, the macro would stop at the first newline.
// Why it is used:
// - it makes the macro behave like one single statement
// - it lets you safely write SUM_FOR(...); with a semicolon after it
// - it avoids strange bugs inside if statements

int main(void) {
    double vec[] = {1.0, 1.0e16, -1.0e16, -0.5};
    size_t n = sizeof(vec) / sizeof(vec[0]);
    double s_for;
    double s_kahan;

    SUM_FOR(s_for, vec, n);         /* Method 1: for loop */
    SUM_KAHAN(s_kahan, vec, n);     /* Method 3: Kahan */

    printf("For loop sum   = %.17g\n", s_for);
#ifdef USE_GSL
    {
        gsl_vector_view v = gsl_vector_view_array(vec, n);
        double s_gsl = gsl_vector_sum(&v.vector);
        printf("GSL sum        = %.17g\n", s_gsl);
    }
#else
    printf("GSL sum        = not compiled (define USE_GSL to enable)\n");
#endif
    printf("Kahan sum      = %.17g\n", s_kahan);
    printf("Exact result   = 0.5\n");
    printf("\n");
    printf("The results may differ because floating-point addition is not associative.\n");
    printf("When 1.0 is added to 1.0e16, the 1.0 can be lost due to rounding.\n");
    printf("For this order, even classic Kahan may still return -0.5.\n");

    return 0;
}
