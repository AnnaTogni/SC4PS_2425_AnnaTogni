// HomeWork6 -- Split the work: daxpy in chunks, partial sums, config file, HDF5
//
//   d[i] = a * x[i] + y[i],   i = 0 ... n-1,   with x[i] = x + i * x_step, y[i] = y + i * y_step
//
//   1. compute d with a single loop (reference) and again chunk by chunk
//      -> check that the two d vectors are identical
//   2. sum every chunk into partial_chunk_sum[], then sum the partial sums
//      -> compare with the sum of d computed by the single loop
//   3. all inputs come from a "Variable = Value" text file (parsed with GLib's GKeyFile)
//   4. all vector chunks and the partial sums are saved in an HDF5 file
//
// ******************************************************
// HOW TO COMPILE:
//   make
// HOW TO RUN:
//   ./chunked_daxpy [config_file]          (default: daxpy.cfg)
//   make test                              (runs the tests in tests/)
// ******************************************************

#define _POSIX_C_SOURCE 200809L // clock_gettime() with -std=c11

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <glib.h>
#include <hdf5.h>

#define GROUP "daxpy" // GKeyFile needs a [group]: we add it in front of the file

typedef struct {
    size_t n;                // vector length
    size_t chunk_size;       // elements per chunk
    size_t number_of_chunks; // ceiling(n / chunk_size)
    double a;
    double x, x_step;        // x[i] = x + i * x_step
    double y, y_step;        // y[i] = y + i * y_step
    char *output;            // HDF5 file name
} Config;

static double seconds_now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

// ****************************************************************************
// 3. Configuration file:  Variable = Value
// ****************************************************************************

// Read a positive integer key. Returns 1 if found, 0 if absent, -1 on error.
static int get_size(GKeyFile *kf, const char *key, size_t *out, const char *path)
{
    if (!g_key_file_has_key(kf, GROUP, key, NULL))
        return 0;
    GError *err = NULL;
    gchar *text = g_key_file_get_value(kf, GROUP, key, &err);
    char *end;
    double v = text ? strtod(text, &end) : 0; // accept 1e6 as well as 1000000
    int ok = text && end != text && *end == '\0' && v >= 1 && v == floor(v) && v <= 1e12;
    if (!ok)
        fprintf(stderr, "%s: '%s' must be a positive integer, got '%s'\n", path, key,
                text ? text : "?");
    else
        *out = (size_t)v;
    g_free(text);
    g_clear_error(&err);
    return ok ? 1 : -1;
}

// Read a double key; 'def' is used if the key is absent and required == 0.
static int get_double(GKeyFile *kf, const char *key, double *out, int required, double def,
                      const char *path)
{
    if (!g_key_file_has_key(kf, GROUP, key, NULL)) {
        if (required) {
            fprintf(stderr, "%s: missing required variable '%s'\n", path, key);
            return 0;
        }
        *out = def;
        return 1;
    }
    GError *err = NULL;
    *out = g_key_file_get_double(kf, GROUP, key, &err);
    if (err != NULL) {
        fprintf(stderr, "%s: '%s' must be a number\n", path, key);
        g_clear_error(&err);
        return 0;
    }
    return 1;
}

static int load_config(const char *path, Config *c)
{
    static const char *known[] = {"n", "chunk_size", "number_of_chunks", "a", "x", "x_step",
                                  "y", "y_step", "output", NULL};
    gchar *raw = NULL;
    GError *err = NULL;
    int ok = 0;

    if (!g_file_get_contents(path, &raw, NULL, &err)) {
        fprintf(stderr, "Cannot read config file: %s\n", err->message);
        g_clear_error(&err);
        return 0;
    }
    gchar *text = g_strdup_printf("[" GROUP "]\n%s", raw);
    GKeyFile *kf = g_key_file_new();
    if (!g_key_file_load_from_data(kf, text, (gsize)-1, G_KEY_FILE_NONE, &err)) {
        fprintf(stderr, "%s: not a valid 'Variable = Value' file: %s\n", path, err->message);
        g_clear_error(&err);
        goto done;
    }

    // Warn about typos: a misspelled variable would otherwise be silently ignored
    gchar **keys = g_key_file_get_keys(kf, GROUP, NULL, NULL);
    for (gchar **k = keys; k && *k; k++) {
        int found = 0;
        for (const char **q = known; *q; q++)
            found |= strcmp(*k, *q) == 0;
        if (!found)
            fprintf(stderr, "%s: warning: unknown variable '%s' ignored\n", path, *k);
    }
    g_strfreev(keys);

    if (get_size(kf, "n", &c->n, path) != 1) {
        if (!g_key_file_has_key(kf, GROUP, "n", NULL))
            fprintf(stderr, "%s: missing required variable 'n'\n", path);
        goto done;
    }

    // Chunks: give either chunk_size or number_of_chunks
    int has_size = get_size(kf, "chunk_size", &c->chunk_size, path);
    size_t nchunks = 0;
    int has_number = get_size(kf, "number_of_chunks", &nchunks, path);
    if (has_size < 0 || has_number < 0)
        goto done;
    if (has_size == has_number) {
        fprintf(stderr, "%s: give exactly one of 'chunk_size' or 'number_of_chunks'\n", path);
        goto done;
    }
    if (has_number) {
        if (nchunks > c->n) {
            fprintf(stderr, "%s: number_of_chunks (%zu) cannot be larger than n (%zu)\n", path,
                    nchunks, c->n);
            goto done;
        }
        c->chunk_size = (c->n + nchunks - 1) / nchunks;
    }
    // CEILING(n / chunk_size) with integers: (n + chunk_size - 1) / chunk_size.
    // Writing ceil(n / chunk_size) in C would NOT work: n / chunk_size between two
    // integers is already truncated before ceil() sees it.
    c->number_of_chunks = (c->n + c->chunk_size - 1) / c->chunk_size;

    if (!get_double(kf, "a", &c->a, 1, 0, path) || !get_double(kf, "x", &c->x, 1, 0, path) ||
        !get_double(kf, "y", &c->y, 1, 0, path) ||
        !get_double(kf, "x_step", &c->x_step, 0, 0.0, path) ||
        !get_double(kf, "y_step", &c->y_step, 0, 0.0, path))
        goto done;

    c->output = g_key_file_get_string(kf, GROUP, "output", NULL);
    if (c->output == NULL || c->output[0] == '\0') {
        g_free(c->output);
        c->output = g_strdup("daxpy_chunks.h5");
    } else if (!g_str_has_suffix(c->output, ".h5") && !g_str_has_suffix(c->output, ".hdf5")) {
        gchar *with_ext = g_strdup_printf("%s.h5", c->output);
        g_free(c->output);
        c->output = with_ext;
    }
    ok = 1;

done:
    g_key_file_free(kf);
    g_free(text);
    g_free(raw);
    return ok;
}

// ****************************************************************************
// 1-2. Computation
// ****************************************************************************

// Reference: a single loop over all elements, and the sum of d in the same loop order.
static double daxpy_reference(const Config *c, const double *x, const double *y, double *d)
{
    double sum = 0.0;
    for (size_t i = 0; i < c->n; i++) {
        d[i] = c->a * x[i] + y[i];
        sum += d[i];
    }
    return sum;
}

// First and one-past-last element of a chunk (0-based: the homework's pseudocode
// is 1-based, current_start = chunk_index * chunk_size + 1, which is off by one in C).
static void chunk_bounds(const Config *c, size_t chunk, size_t *start, size_t *end)
{
    *start = chunk * c->chunk_size;
    *end = *start + c->chunk_size < c->n ? *start + c->chunk_size : c->n;
}

// Chunked: outer loop over chunks, inner loop over the elements of one chunk.
static double daxpy_chunked(const Config *c, const double *x, const double *y, double *d,
                            double *partial_chunk_sum)
{
    for (size_t chunk = 0; chunk < c->number_of_chunks; chunk++) {
        size_t start, end;
        chunk_bounds(c, chunk, &start, &end);
        double partial = 0.0;
        for (size_t i = start; i < end; i++) {
            d[i] = c->a * x[i] + y[i];
            partial += d[i];
        }
        partial_chunk_sum[chunk] = partial;
    }
    double total = 0.0;
    for (size_t chunk = 0; chunk < c->number_of_chunks; chunk++)
        total += partial_chunk_sum[chunk];
    return total;
}

// Neumaier compensated sum (see HomeWork5): a much more accurate estimate of the sum,
// used to judge which of the two sums above is closer to the truth.
static double sum_neumaier(const double *v, size_t n)
{
    double s = 0.0, c = 0.0;
    for (size_t i = 0; i < n; i++) {
        double t = s + v[i];
        c += fabs(s) >= fabs(v[i]) ? (s - t) + v[i] : (v[i] - t) + s;
        s = t;
    }
    return s + c;
}

// ****************************************************************************
// 4. HDF5 output
//
//   /                      attributes: n, chunk_size, number_of_chunks, a
//   /x, /y, /d             full vectors (length n)
//   /partial_chunk_sum     one value per chunk
//   /total_sum             attributes: chunked, reference, accurate (Neumaier)
//   /chunks/x, /y, /d      2-D datasets [number_of_chunks][chunk_size]: row k = chunk k
//                          (the last row is padded with NaN if the last chunk is shorter)
//   /chunks/start, /length first element and number of elements of every chunk
//
// One row per chunk (instead of one group per chunk) keeps the file small and fast
// also with millions of chunks.
// ****************************************************************************

static int write_vector(hid_t where, const char *name, const double *data, size_t len)
{
    hsize_t dims[1] = {len};
    hid_t space = H5Screate_simple(1, dims, NULL);
    hid_t set = space >= 0 ? H5Dcreate2(where, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT,
                                        H5P_DEFAULT, H5P_DEFAULT)
                           : -1;
    herr_t status = set >= 0 ? H5Dwrite(set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)
                             : -1;
    if (set >= 0)
        H5Dclose(set);
    if (space >= 0)
        H5Sclose(space);
    return status >= 0;
}

static int write_index(hid_t where, const char *name, const unsigned long long *data, size_t len)
{
    hsize_t dims[1] = {len};
    hid_t space = H5Screate_simple(1, dims, NULL);
    hid_t set = space >= 0 ? H5Dcreate2(where, name, H5T_NATIVE_ULLONG, space, H5P_DEFAULT,
                                        H5P_DEFAULT, H5P_DEFAULT)
                           : -1;
    herr_t status = set >= 0 ? H5Dwrite(set, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, data)
                             : -1;
    if (set >= 0)
        H5Dclose(set);
    if (space >= 0)
        H5Sclose(space);
    return status >= 0;
}

// Write a vector of length n as a 2-D dataset [number_of_chunks][chunk_size], row k = chunk k.
// The vector is contiguous and every full chunk is a full row, so the full rows are written
// with one H5Dwrite; the (possibly shorter) last chunk is written separately.
static int write_chunked(hid_t where, const char *name, const double *v, const Config *c)
{
    hsize_t dims[2] = {c->number_of_chunks, c->chunk_size};
    double nan = NAN;
    hid_t space = H5Screate_simple(2, dims, NULL);
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    int ok = space >= 0 && plist >= 0 && H5Pset_fill_value(plist, H5T_NATIVE_DOUBLE, &nan) >= 0;
    hid_t set = ok ? H5Dcreate2(where, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, plist,
                                H5P_DEFAULT)
                   : -1;
    ok = set >= 0;

    size_t last_len = c->n - (c->number_of_chunks - 1) * c->chunk_size;
    size_t full_rows = last_len == c->chunk_size ? c->number_of_chunks : c->number_of_chunks - 1;
    hsize_t offset[2] = {0, 0}, count[2] = {full_rows, c->chunk_size};
    if (ok && full_rows > 0) {
        hid_t mem = H5Screate_simple(2, count, NULL);
        ok = mem >= 0 &&
             H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, count, NULL) >= 0 &&
             H5Dwrite(set, H5T_NATIVE_DOUBLE, mem, space, H5P_DEFAULT, v) >= 0;
        if (mem >= 0)
            H5Sclose(mem);
    }
    if (ok && full_rows < c->number_of_chunks) {
        hsize_t off2[2] = {c->number_of_chunks - 1, 0}, cnt2[2] = {1, last_len};
        hid_t mem = H5Screate_simple(2, cnt2, NULL);
        ok = mem >= 0 && H5Sselect_hyperslab(space, H5S_SELECT_SET, off2, NULL, cnt2, NULL) >= 0 &&
             H5Dwrite(set, H5T_NATIVE_DOUBLE, mem, space, H5P_DEFAULT,
                      v + (c->number_of_chunks - 1) * c->chunk_size) >= 0;
        if (mem >= 0)
            H5Sclose(mem);
    }
    if (set >= 0)
        H5Dclose(set);
    if (plist >= 0)
        H5Pclose(plist);
    if (space >= 0)
        H5Sclose(space);
    return ok;
}

static int write_attr(hid_t where, const char *name, hid_t type, const void *value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = space >= 0 ? H5Acreate2(where, name, type, space, H5P_DEFAULT, H5P_DEFAULT) : -1;
    herr_t status = attr >= 0 ? H5Awrite(attr, type, value) : -1;
    if (attr >= 0)
        H5Aclose(attr);
    if (space >= 0)
        H5Sclose(space);
    return status >= 0;
}

static int write_hdf5(const Config *c, const double *x, const double *y, const double *d,
                      const double *partial, double total_chunked, double total_reference,
                      double total_accurate)
{
    hid_t file = H5Fcreate(c->output, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0)
        return 0;

    unsigned long long n = c->n, cs = c->chunk_size, nc = c->number_of_chunks;
    int ok = write_attr(file, "n", H5T_NATIVE_ULLONG, &n) &&
             write_attr(file, "chunk_size", H5T_NATIVE_ULLONG, &cs) &&
             write_attr(file, "number_of_chunks", H5T_NATIVE_ULLONG, &nc) &&
             write_attr(file, "a", H5T_NATIVE_DOUBLE, &c->a) && write_vector(file, "x", x, c->n) &&
             write_vector(file, "y", y, c->n) && write_vector(file, "d", d, c->n) &&
             write_vector(file, "partial_chunk_sum", partial, c->number_of_chunks);

    if (ok) {
        hid_t g = H5Gcreate2(file, "total_sum", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ok = g >= 0 && write_attr(g, "chunked", H5T_NATIVE_DOUBLE, &total_chunked) &&
             write_attr(g, "reference", H5T_NATIVE_DOUBLE, &total_reference) &&
             write_attr(g, "accurate", H5T_NATIVE_DOUBLE, &total_accurate);
        if (g >= 0)
            H5Gclose(g);
    }

    if (ok) {
        unsigned long long *start = malloc(c->number_of_chunks * sizeof(*start));
        unsigned long long *length = malloc(c->number_of_chunks * sizeof(*length));
        hid_t g = H5Gcreate2(file, "chunks", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ok = start && length && g >= 0;
        for (size_t k = 0; ok && k < c->number_of_chunks; k++) {
            size_t s, e;
            chunk_bounds(c, k, &s, &e);
            start[k] = s;
            length[k] = e - s;
        }
        ok = ok && write_chunked(g, "x", x, c) && write_chunked(g, "y", y, c) &&
             write_chunked(g, "d", d, c) && write_index(g, "start", start, c->number_of_chunks) &&
             write_index(g, "length", length, c->number_of_chunks);
        if (g >= 0)
            H5Gclose(g);
        free(start);
        free(length);
    }
    return H5Fclose(file) >= 0 && ok;
}

// ****************************************************************************

int main(int argc, char *argv[])
{
    if (argc > 2 || (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))) {
        fprintf(stderr, "Usage: %s [config_file]   (default: daxpy.cfg)\n", argv[0]);
        return 1;
    }
    const char *path = argc == 2 ? argv[1] : "daxpy.cfg";

    Config c = {0};
    if (!load_config(path, &c))
        return 1;

    double *x = malloc(c.n * sizeof(*x));
    double *y = malloc(c.n * sizeof(*y));
    double *d_ref = malloc(c.n * sizeof(*d_ref));
    double *d_chunk = malloc(c.n * sizeof(*d_chunk));
    double *partial = malloc(c.number_of_chunks * sizeof(*partial));
    int status = 1;
    if (!x || !y || !d_ref || !d_chunk || !partial) {
        fprintf(stderr, "Error: not enough memory for n = %zu\n", c.n);
        goto cleanup;
    }
    for (size_t i = 0; i < c.n; i++) {
        x[i] = c.x + (double)i * c.x_step;
        y[i] = c.y + (double)i * c.y_step;
    }

    printf("Config file:      %s\n", path);
    printf("n = %zu, chunk_size = %zu, number_of_chunks = %zu", c.n, c.chunk_size,
           c.number_of_chunks);
    if (c.n % c.chunk_size != 0)
        printf(" (last chunk: %zu elements)", c.n - (c.number_of_chunks - 1) * c.chunk_size);
    printf("\na = %g, x[i] = %g + i*%g, y[i] = %g + i*%g\n\n", c.a, c.x, c.x_step, c.y, c.y_step);

    // ---- 1. single loop vs chunks --------------------------------------------
    double t0 = seconds_now();
    double sum_ref = daxpy_reference(&c, x, y, d_ref);
    double t1 = seconds_now();
    double sum_chunk = daxpy_chunked(&c, x, y, d_chunk, partial);
    double t2 = seconds_now();

    size_t n_diff = 0;
    for (size_t i = 0; i < c.n; i++)
        n_diff += d_ref[i] != d_chunk[i];
    printf("1. d (chunks) == d (single loop), element by element: %s (%zu of %zu differ)\n",
           n_diff == 0 ? "YES" : "NO", n_diff, c.n);
    printf("   time: single loop %.3e s, chunked %.3e s\n\n", t1 - t0, t2 - t1);

    // ---- 2. partial sums --------------------------------------------------------
    // The two sums add the SAME numbers in a DIFFERENT order, so they can differ by
    // rounding. Error bound of a sum of n terms: about n * eps * sum|d_i|.
    double sum_abs = 0.0;
    for (size_t i = 0; i < c.n; i++)
        sum_abs += fabs(d_ref[i]);
    double tol = (double)c.n * DBL_EPSILON * sum_abs;
    double accurate = sum_neumaier(d_ref, c.n);
    int sums_ok = fabs(sum_chunk - sum_ref) <= tol;

    printf("2. sum of d\n");
    printf("   single loop:                 %.17g   (error %.2e)\n", sum_ref, sum_ref - accurate);
    printf("   sum of partial_chunk_sum[]:  %.17g   (error %.2e)\n", sum_chunk, sum_chunk - accurate);
    printf("   accurate (Neumaier):         %.17g\n", accurate);
    printf("   exactly equal: %s;  equal within rounding (|diff| <= n eps sum|d| = %.1e): %s\n\n",
           sum_chunk == sum_ref ? "YES" : "NO", tol, sums_ok ? "YES" : "NO");

    // ---- 4. HDF5 ------------------------------------------------------------------
    t0 = seconds_now();
    int written = write_hdf5(&c, x, y, d_chunk, partial, sum_chunk, sum_ref, accurate);
    t1 = seconds_now();
    if (written)
        printf("4. HDF5 file '%s' written in %.3e s (inspect: h5ls -r %s)\n", c.output, t1 - t0,
               c.output);
    else
        fprintf(stderr, "4. Error: could not write the HDF5 file '%s'\n", c.output);

    status = (n_diff == 0 && sums_ok && written) ? 0 : 2;

cleanup:
    free(x);
    free(y);
    free(d_ref);
    free(d_chunk);
    free(partial);
    g_free(c.output);
    return status;
}
