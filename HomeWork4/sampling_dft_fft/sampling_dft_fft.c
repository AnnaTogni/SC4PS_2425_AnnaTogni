// HomeWork4 -- Sampling, DFT, FFT and aliasing
//
// C version of sampling_dft_fft_alias_macro_old.C (ROOT / C++).
//
// A continuous signal  x(t) = sum_i sin(2 pi f_i t)  is sampled at frequency fs
// for a duration T (N = fs T samples). The program then:
//   1. predicts where every input frequency is observed after sampling (aliasing)
//   2. computes the spectrum with a naive DFT, O(N^2), and with an FFT, O(N log N)
//   3. reconstructs the samples with the inverse transforms
//   4. times both methods and checks that they agree
//   5. writes CSV files that plot_sampling_dft_fft.gp turns into a figure
//
// The FFT is a recursive mixed-radix Cooley-Tukey algorithm: it works for ANY N
// (not only powers of two), so the requested duration T is sampled exactly.
//
// ******************************************************
// HOW TO COMPILE:
//   make
// HOW TO RUN:
//   ./sampling_dft_fft fs T f1 [f2 ...]
//   ./sampling_dft_fft 128 1 50 120       (120 Hz aliases to 8 Hz)
//   make plot                             (figure in plots/sampling_dft_fft.png)
// ******************************************************

#define _POSIX_C_SOURCE 200809L // for clock_gettime() and mkdir() with -std=c11

#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define OUTPUT_DIR   "output"
#define MAX_SAMPLES  (1u << 22) // ~4 million samples (64 MB per complex array)
#define MAX_DFT      16384u     // above this the O(N^2) DFT would take too long
#define PEAK_REL     0.05       // a spectral peak must be > 5% of the largest one

typedef double complex cplx;

// ****************************************************************************
// Input helpers
// ****************************************************************************

static void usage(const char *prog)
{
    fprintf(stderr,
            "Usage: %s fs T f1 [f2 ...]\n"
            "  fs   sampling frequency [Hz]\n"
            "  T    duration of the sampled signal [s]  (N = fs * T samples)\n"
            "  f_i  frequencies [Hz] of the sines in x(t) = sum_i sin(2 pi f_i t)\n"
            "Example: %s 128 1 50 120\n",
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

static int ensure_directory(const char *path)
{
    if (mkdir(path, 0775) == 0 || errno == EEXIST)
        return 0;
    perror(path);
    return -1;
}

// ****************************************************************************
// Signal and aliasing
// ****************************************************************************

static double signal_value(double t, const double *freqs, size_t n_freqs)
{
    double value = 0.0;
    for (size_t i = 0; i < n_freqs; i++)
        value += sin(2.0 * M_PI * freqs[i] * t);
    return value;
}

// Frequency at which f is observed after sampling at fs: fold f into [0, fs/2].
static double alias_frequency(double f, double fs)
{
    double wrapped = fmod(fabs(f), fs);
    if (wrapped > fs / 2.0)
        wrapped = fs - wrapped;
    return wrapped;
}

// ****************************************************************************
// Transforms. Convention:
//   forward  X[k] = sum_m x[m] exp(-2 pi i k m / N)
//   inverse  x[m] = 1/N sum_k X[k] exp(+2 pi i k m / N)
// Both methods use the same table of "twiddle factors" w[j] = exp(sign 2 pi i j / N),
// so the timing compares the algorithms, not the cost of computing sin/cos.
// ****************************************************************************

static cplx *twiddle_table(size_t n, int sign)
{
    cplx *w = malloc(n * sizeof(*w));
    if (w == NULL)
        return NULL;
    for (size_t j = 0; j < n; j++)
        w[j] = cexp(sign * 2.0 * M_PI * I * (double)j / (double)n);
    return w;
}

// Naive DFT: N^2 multiplications.
static void dft(const cplx *in, cplx *out, size_t n, const cplx *w)
{
    for (size_t k = 0; k < n; k++) {
        cplx sum = 0.0;
        for (size_t m = 0; m < n; m++)
            sum += in[m] * w[(k * m) % n];
        out[k] = sum;
    }
}

static size_t smallest_factor(size_t n)
{
    if (n % 2 == 0)
        return 2;
    for (size_t p = 3; p * p <= n; p += 2)
        if (n % p == 0)
            return p;
    return n; // n is prime
}

// Mixed-radix Cooley-Tukey FFT of the n values in[0], in[stride], in[2 stride], ...
// n_total is the length of the full transform (the twiddle table has n_total entries);
// tmp is a work array of at least n elements.
//
// If n = p * m (p = smallest prime factor), the n inputs are split into p
// interleaved sub-sequences of length m, each is transformed recursively, and
// the p results are combined:
//   X[k + q m] = sum_{r=0}^{p-1} W_n^{r (k + q m)} Y_r[k],   W_n = exp(sign 2 pi i / n)
// For n = 2^k this is the classic radix-2 FFT. For a prime n it is a direct DFT.
static void fft_rec(const cplx *in, size_t stride, cplx *out, size_t n,
                    const cplx *w, size_t n_total, cplx *tmp)
{
    if (n == 1) {
        out[0] = in[0];
        return;
    }

    size_t p = smallest_factor(n);
    size_t m = n / p;
    size_t w_step = n_total / n; // W_n^j = w[j * w_step]

    if (p == n) { // prime length: direct DFT
        for (size_t k = 0; k < n; k++) {
            cplx sum = 0.0;
            for (size_t j = 0; j < n; j++)
                sum += in[j * stride] * w[((j * k) % n) * w_step];
            out[k] = sum;
        }
        return;
    }

    // Sub-transform r goes into out[r m ... r m + m - 1].
    for (size_t r = 0; r < p; r++)
        fft_rec(in + r * stride, stride * p, out + r * m, m, w, n_total, tmp);

    for (size_t q = 0; q < p; q++)
        for (size_t k = 0; k < m; k++) {
            size_t kk = k + q * m;
            cplx sum = 0.0;
            for (size_t r = 0; r < p; r++)
                sum += out[r * m + k] * w[((r * kk) % n) * w_step];
            tmp[kk] = sum;
        }
    memcpy(out, tmp, n * sizeof(*out));
}

// ****************************************************************************
// Spectrum analysis
// ****************************************************************************

// One-sided amplitude spectrum: a sine of amplitude A on a bin gives a peak A.
static void amplitude_spectrum(const cplx *X, size_t n, double *amp)
{
    for (size_t k = 0; k <= n / 2; k++) {
        double scale = (k == 0 || (n % 2 == 0 && k == n / 2)) ? 1.0 : 2.0;
        amp[k] = scale * cabs(X[k]) / (double)n;
    }
}

static double max_abs_diff_c(const cplx *a, const cplx *b, size_t n)
{
    double d = 0.0;
    for (size_t i = 0; i < n; i++)
        d = fmax(d, cabs(a[i] - b[i]));
    return d;
}

static double max_abs_diff_r(const double *a, const cplx *b, size_t n)
{
    double d = 0.0;
    for (size_t i = 0; i < n; i++)
        d = fmax(d, cabs(a[i] - b[i]));
    return d;
}

// Print the local maxima of the spectrum that are > PEAK_REL of the largest one.
static void report_peaks(const double *amp, size_t bins, double df)
{
    double max_amp = 0.0;
    for (size_t k = 0; k < bins; k++)
        max_amp = fmax(max_amp, amp[k]);

    if (max_amp < 1e-9) {
        printf("  none: the spectrum is zero (all the samples are zero)\n");
        return;
    }
    for (size_t k = 0; k < bins; k++) {
        double left = k > 0 ? amp[k - 1] : 0.0;
        double right = k + 1 < bins ? amp[k + 1] : 0.0;
        if (amp[k] > PEAK_REL * max_amp && amp[k] >= left && amp[k] > right)
            printf("  f = %10.4f Hz   amplitude = %.4f\n", k * df, amp[k]);
    }
}

// ****************************************************************************

int main(int argc, char *argv[])
{
    // ---- 1. Input ----------------------------------------------------------
    if (argc < 4) {
        usage(argv[0]);
        return 1;
    }
    double fs, duration;
    if (!parse_double(argv[1], &fs) || fs <= 0.0 ||
        !parse_double(argv[2], &duration) || duration <= 0.0) {
        fprintf(stderr, "Error: fs and T must be positive numbers\n");
        usage(argv[0]);
        return 1;
    }
    size_t n_freqs = (size_t)(argc - 3);
    double *freqs = malloc(n_freqs * sizeof(*freqs));
    if (freqs == NULL) {
        fprintf(stderr, "Error: not enough memory\n");
        return 1;
    }
    for (size_t i = 0; i < n_freqs; i++)
        if (!parse_double(argv[3 + i], &freqs[i]) || freqs[i] < 0.0) {
            fprintf(stderr, "Error: frequencies must be numbers >= 0, got '%s'\n", argv[3 + i]);
            free(freqs);
            return 1;
        }

    double n_real = round(fs * duration);
    if (n_real < 8 || n_real > MAX_SAMPLES) {
        fprintf(stderr, "Error: N = fs * T = %.0f samples, must be between 8 and %u\n",
                n_real, MAX_SAMPLES);
        free(freqs);
        return 1;
    }
    size_t N = (size_t)n_real;
    size_t bins = N / 2 + 1;
    double dt = 1.0 / fs;
    double f_nyquist = fs / 2.0;
    double df = fs / (double)N;

    printf("===== Sampling =====\n");
    printf("fs = %g Hz, dt = %g s, N = %zu samples, T = N dt = %g s\n", fs, dt, N, N * dt);
    printf("Nyquist frequency  f_N = fs/2 = %g Hz\n", f_nyquist);
    printf("Frequency step     df  = fs/N = 1/T = %g Hz\n", df);
    if (fabs(N * dt - duration) > 1e-9 * duration)
        printf("(fs * T is not an integer: T rounded to %g s)\n", N * dt);

    // ---- 2. Where is every input frequency observed? -----------------------
    printf("\n===== Input frequencies and aliasing =====\n");
    for (size_t i = 0; i < n_freqs; i++) {
        double f = freqs[i];
        double fa = alias_frequency(f, fs);
        printf("f = %10.4f Hz -> ", f);
        // A sine whose observed frequency is exactly 0 or fs/2 is sampled at its
        // zeros: sin(2 pi fa n / fs) = sin(0) or sin(pi n) = 0 for every n.
        if (fa < 1e-9 * fs || fabs(fa - f_nyquist) < 1e-9 * fs)
            printf("observed at %g Hz, where a sine is sampled at its zeros:"
                   " the component is INVISIBLE\n", fa);
        else if (f < f_nyquist)
            printf("below Nyquist, observed at its true frequency\n");
        else
            printf("ABOVE Nyquist, aliases to %.4f Hz\n", fa);
        if (fabs(fa / df - round(fa / df)) > 1e-6 && fa > 0 && fa < f_nyquist)
            printf("%17s(not a multiple of df = %g Hz: expect leakage into nearby bins)\n",
                   "", df);
    }

    // ---- 3. Sample the signal ------------------------------------------------
    double *x = malloc(N * sizeof(*x));
    cplx *xc = malloc(N * sizeof(*xc));
    cplx *X_dft = malloc(N * sizeof(*X_dft));
    cplx *X_fft = malloc(N * sizeof(*X_fft));
    cplx *x_dft = malloc(N * sizeof(*x_dft));
    cplx *x_fft = malloc(N * sizeof(*x_fft));
    cplx *tmp = malloc(N * sizeof(*tmp));
    double *amp_dft = malloc(bins * sizeof(*amp_dft));
    double *amp_fft = malloc(bins * sizeof(*amp_fft));
    cplx *w_fwd = twiddle_table(N, -1);
    cplx *w_inv = twiddle_table(N, +1);
    int status = 1;
    if (!x || !xc || !X_dft || !X_fft || !x_dft || !x_fft || !tmp || !amp_dft ||
        !amp_fft || !w_fwd || !w_inv) {
        fprintf(stderr, "Error: not enough memory for N = %zu\n", N);
        goto cleanup;
    }

    for (size_t m = 0; m < N; m++) {
        x[m] = signal_value(m * dt, freqs, n_freqs);
        xc[m] = x[m];
    }

    // ---- 4. Transforms and timing -----------------------------------------
    int do_dft = N <= MAX_DFT;
    double t_dft = 0.0, t_idft = 0.0;
    if (do_dft) {
        double t0 = seconds_now();
        dft(xc, X_dft, N, w_fwd);
        double t1 = seconds_now();
        dft(X_dft, x_dft, N, w_inv);
        double t2 = seconds_now();
        for (size_t m = 0; m < N; m++)
            x_dft[m] /= (double)N;
        t_dft = t1 - t0;
        t_idft = t2 - t1;
    }

    double t0 = seconds_now();
    fft_rec(xc, 1, X_fft, N, w_fwd, N, tmp);
    double t1 = seconds_now();
    fft_rec(X_fft, 1, x_fft, N, w_inv, N, tmp);
    double t2 = seconds_now();
    for (size_t m = 0; m < N; m++)
        x_fft[m] /= (double)N;
    double t_fft = t1 - t0, t_ifft = t2 - t1;

    printf("\n===== Timing (N = %zu) =====\n", N);
    if (do_dft) {
        printf("naive DFT: forward %10.3f ms   inverse %10.3f ms\n", 1e3 * t_dft, 1e3 * t_idft);
        printf("FFT:       forward %10.3f ms   inverse %10.3f ms   -> %.0fx faster\n",
               1e3 * t_fft, 1e3 * t_ifft, (t_dft + t_idft) / (t_fft + t_ifft));
    } else {
        printf("naive DFT: skipped (N > %u would take too long: it scales as N^2)\n", MAX_DFT);
        printf("FFT:       forward %10.3f ms   inverse %10.3f ms\n", 1e3 * t_fft, 1e3 * t_ifft);
    }
    printf("prime factor used first by the FFT: %zu%s\n", smallest_factor(N),
           smallest_factor(N) == N ? " (N is prime: the FFT is a plain DFT!)" : "");

    // ---- 5. Checks -----------------------------------------------------------
    double x_max = 0.0;
    for (size_t m = 0; m < N; m++)
        x_max = fmax(x_max, fabs(x[m]));
    double tol = 1e-9 * fmax(1.0, x_max) * N; // rounding grows with N

    printf("\n===== Checks =====\n");
    double err_ifft = max_abs_diff_r(x, x_fft, N);
    int ok = err_ifft <= tol;
    printf("max |x - IFFT(FFT(x))|  = %.2e   %s\n", err_ifft, err_ifft <= tol ? "OK" : "FAILED");
    if (do_dft) {
        double err_idft = max_abs_diff_r(x, x_dft, N);
        double err_spec = max_abs_diff_c(X_dft, X_fft, N);
        printf("max |x - IDFT(DFT(x))|  = %.2e   %s\n", err_idft, err_idft <= tol ? "OK" : "FAILED");
        printf("max |DFT(x) - FFT(x)|   = %.2e   %s\n", err_spec, err_spec <= tol * N ? "OK" : "FAILED");
        ok = ok && err_idft <= tol && err_spec <= tol * N;
    }

    amplitude_spectrum(X_fft, N, amp_fft);
    if (do_dft)
        amplitude_spectrum(X_dft, N, amp_dft);

    printf("\n===== Peaks in the FFT spectrum (0 ... %g Hz) =====\n", f_nyquist);
    report_peaks(amp_fft, bins, df);

    // ---- 6. Output files for the plots -------------------------------------
    if (ensure_directory(OUTPUT_DIR) != 0)
        goto cleanup;

    FILE *f;
    // Continuous signal on a fine grid (20 points per sample, at most 200000).
    size_t n_fine = N * 20 < 200000 ? N * 20 : 200000;
    if ((f = fopen(OUTPUT_DIR "/continuous.csv", "w")) == NULL)
        goto file_error;
    fprintf(f, "time_s,signal\n");
    for (size_t i = 0; i <= n_fine; i++) {
        double t = N * dt * (double)i / (double)n_fine;
        fprintf(f, "%.10g,%.10g\n", t, signal_value(t, freqs, n_freqs));
    }
    fclose(f);

    if ((f = fopen(OUTPUT_DIR "/samples.csv", "w")) == NULL)
        goto file_error;
    fprintf(f, "time_s,sample,reconstructed_dft,reconstructed_fft\n");
    for (size_t m = 0; m < N; m++)
        fprintf(f, "%.10g,%.10g,%.10g,%.10g\n", m * dt, x[m],
                do_dft ? creal(x_dft[m]) : NAN, creal(x_fft[m]));
    fclose(f);

    if ((f = fopen(OUTPUT_DIR "/spectrum.csv", "w")) == NULL)
        goto file_error;
    fprintf(f, "frequency_hz,amplitude_dft,amplitude_fft\n");
    for (size_t k = 0; k < bins; k++)
        fprintf(f, "%.10g,%.10g,%.10g\n", k * df, do_dft ? amp_dft[k] : NAN, amp_fft[k]);
    fclose(f);

    if ((f = fopen(OUTPUT_DIR "/frequencies.csv", "w")) == NULL)
        goto file_error;
    fprintf(f, "input_hz,observed_hz,nyquist_hz\n");
    for (size_t i = 0; i < n_freqs; i++)
        fprintf(f, "%.10g,%.10g,%.10g\n", freqs[i], alias_frequency(freqs[i], fs), f_nyquist);
    fclose(f);

    printf("\nData written to %s/ (continuous, samples, spectrum, frequencies .csv)\n",
           OUTPUT_DIR);
    status = ok ? 0 : 2;
    goto cleanup;

file_error:
    perror("Error writing the output files");

cleanup:
    free(freqs);
    free(x);
    free(xc);
    free(X_dft);
    free(X_fft);
    free(x_dft);
    free(x_fft);
    free(tmp);
    free(amp_dft);
    free(amp_fft);
    free(w_fwd);
    free(w_inv);
    return status;
}
