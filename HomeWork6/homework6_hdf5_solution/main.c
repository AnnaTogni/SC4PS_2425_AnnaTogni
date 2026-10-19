#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <glib.h>

// Check for HDF5 headers in common locations. 
// We want to support both the standard hdf5.h and the one installed by some package managers under hdf5/serial/hdf5.h.
#if defined(__has_include)
#if __has_include(<hdf5.h>)
#include <hdf5.h>
#define HAVE_HDF5 1
#elif __has_include(<hdf5/serial/hdf5.h>)
#include <hdf5/serial/hdf5.h>
#define HAVE_HDF5 1
#endif
#endif

// If we didn't find HDF5 headers, define HAVE_HDF5 as 0 to disable HDF5-related code.
#ifndef HAVE_HDF5
#define HAVE_HDF5 0
#endif

// Configuration structure to hold parameters read from the config file or user input.
typedef struct {
    gsize vector_length;
    gsize chunk_size;
    double alpha;
    gchar *output_hdf5;
} Config;

// DataSet structure to hold the vectors and intermediate results.
typedef struct {
    double *x;
    double *y;
    double *d_reference;
    double *d_chunked;
    double *partial_chunk_sum;
    gsize number_of_chunks;
} DataSet;

static gchar *normalize_hdf5_filename(const char *input_name);
static double wall_time_seconds(void);

// Helper functions for memory management, file reading, configuration loading, vector initialization, computation, comparison, and HDF5 output.
static void free_config(Config *config) {
    g_free(config->output_hdf5);
}

// Frees all dynamically allocated arrays in the DataSet structure.
static void free_dataset(DataSet *data) {
    g_free(data->x);
    g_free(data->y);
    g_free(data->d_reference);
    g_free(data->d_chunked);
    g_free(data->partial_chunk_sum);
}

// Reads the entire contents of a text file into a string. 
// The caller is responsible for freeing the content string.
static gboolean read_text_file(const gchar *path, gchar **content, GError **error) {
    return g_file_get_contents(path, content, NULL, error);
}

// Loads the configuration parameters from a file. 
// The file is expected to have key-value pairs that can be parsed by GKeyFile.
static gboolean load_config(const gchar *path, Config *config, GError **error) {
    gchar *raw_content = NULL;
    gchar *wrapped_content = NULL;
    GKeyFile *key_file = NULL;
    gboolean ok = FALSE;

    if (!read_text_file(path, &raw_content, error)) {
        return FALSE;
    }

    wrapped_content = g_strdup_printf("[homework]\n%s", raw_content);
    key_file = g_key_file_new();

    if (!g_key_file_load_from_data(key_file, wrapped_content, -1, G_KEY_FILE_NONE, error)) {
        goto cleanup;
    }

    config->vector_length = (gsize)g_key_file_get_uint64(key_file, "homework", "vector_length", error);
    if (*error != NULL) {
        goto cleanup;
    }

    config->chunk_size = (gsize)g_key_file_get_uint64(key_file, "homework", "chunk_size", error);
    if (*error != NULL) {
        goto cleanup;
    }

    config->alpha = g_key_file_get_double(key_file, "homework", "alpha", error);
    if (*error != NULL) {
        goto cleanup;
    }

    config->output_hdf5 = g_key_file_get_string(key_file, "homework", "output_hdf5", error);
    if (*error != NULL) {
        goto cleanup;
    }

    {
        gchar *normalized_name = normalize_hdf5_filename(config->output_hdf5);
        g_free(config->output_hdf5);
        config->output_hdf5 = normalized_name;
    }

    if (config->vector_length == 0 || config->chunk_size == 0) {
        g_set_error(error, G_KEY_FILE_ERROR, G_KEY_FILE_ERROR_INVALID_VALUE,
                    "vector_length and chunk_size must be greater than zero");
        goto cleanup;
    }

    ok = TRUE;

cleanup:
    g_key_file_unref(key_file);
    g_free(wrapped_content);
    g_free(raw_content);
    return ok;
}

// Allocates memory for the vectors and intermediate results based on the configuration parameters.
static gboolean allocate_dataset(const Config *config, DataSet *data, GError **error) {
    data->number_of_chunks = (config->vector_length + config->chunk_size - 1) / config->chunk_size;

    data->x = g_new0(double, config->vector_length);
    data->y = g_new0(double, config->vector_length);
    data->d_reference = g_new0(double, config->vector_length);
    data->d_chunked = g_new0(double, config->vector_length);
    data->partial_chunk_sum = g_new0(double, data->number_of_chunks);

    if (data->x == NULL || data->y == NULL || data->d_reference == NULL ||
        data->d_chunked == NULL || data->partial_chunk_sum == NULL) {
        g_set_error(error, G_FILE_ERROR, g_file_error_from_errno(ENOMEM),
                    "failed to allocate memory");
        return FALSE;
    }

    return TRUE;
}

// Helper functions to read different types of values from standard input with prompts.
static gboolean read_gsize_value(const char *prompt, gsize *value) {
    size_t temp = 0;
    printf("%s", prompt);
    if (scanf("%zu", &temp) != 1 || temp == 0) {
        return FALSE;
    }
    *value = (gsize)temp;
    return TRUE;
}

// Reads a double value from standard input with a prompt. 
// Returns TRUE on success, FALSE on failure.
static gboolean read_double_value(const char *prompt, double *value) {
    printf("%s", prompt);
    return scanf("%lf", value) == 1;
}

// Reads a character value from standard input with a prompt. 
// Returns TRUE on success, FALSE on failure.
static gboolean read_char_value(const char *prompt, char *value) {
    printf("%s", prompt);
    return scanf(" %c", value) == 1;
}

// Reads a string value from standard input with a prompt. 
// Returns TRUE on success, FALSE on failure.
static gboolean read_string_value(const char *prompt, char *buffer, size_t buffer_size) {
    printf("%s", prompt);
    return scanf("%255s", buffer) == 1 && strlen(buffer) < buffer_size;
}

static gchar *normalize_hdf5_filename(const char *input_name) {
    if (g_str_has_suffix(input_name, ".h5") || g_str_has_suffix(input_name, ".hdf5")) {
        return g_strdup(input_name);
    }

    return g_strdup_printf("%s.h5", input_name);
}

static double wall_time_seconds(void) {
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1.0e9;
}

// Fills a vector interactively by asking the user if they want all elements to be the same or to enter them individually.
static gboolean fill_vector_interactively(double *vector, gsize length, const char *name) {
    char answer;
    double value;

    printf("Do you want all elements of %s to be the same? (y/n): ", name);
    if (!read_char_value("", &answer)) {
        return FALSE;
    }

    if (answer == 'y' || answer == 'Y') {
        char prompt[64];
        snprintf(prompt, sizeof(prompt), "Value for %s: ", name);
        if (!read_double_value(prompt, &value)) {
            return FALSE;
        }

        for (gsize i = 0; i < length; ++i) {
            vector[i] = value;
        }
    } else {
        for (gsize i = 0; i < length; ++i) {
            printf("%s[%zu] = ", name, i);
            if (scanf("%lf", &vector[i]) != 1) {
                return FALSE;
            }
        }
    }

    return TRUE;
}

// Initializes the vectors x and y by asking the user for input.
static gboolean initialize_vectors_from_input(const Config *config, DataSet *data) {
    if (!fill_vector_interactively(data->x, config->vector_length, "x")) {
        return FALSE;
    }

    if (!fill_vector_interactively(data->y, config->vector_length, "y")) {
        return FALSE;
    }

    return TRUE;
}

// Initializes the vectors x and y by reading values from the configuration file using GKeyFile.
static void initialize_vectors_from_config(const Config *config, DataSet *data, GKeyFile *key_file) {
    for (gsize i = 0; i < config->vector_length; ++i) {
        gchar *x_key = g_strdup_printf("x_%zu", i);
        gchar *y_key = g_strdup_printf("y_%zu", i);
        data->x[i] = g_key_file_get_double(key_file, "homework", x_key, NULL);
        data->y[i] = g_key_file_get_double(key_file, "homework", y_key, NULL);
        g_free(x_key);
        g_free(y_key);
    }
}

// Loads the configuration parameters and initializes the vectors from a configuration file.
static gboolean load_config_and_vectors(const gchar *path, Config *config, DataSet *data, GError **error) {
    gchar *raw_content = NULL;
    gchar *wrapped_content = NULL;
    GKeyFile *key_file = NULL;
    gboolean ok = FALSE;

    if (!load_config(path, config, error)) {
        return FALSE;
    }

    if (!allocate_dataset(config, data, error)) {
        return FALSE;
    }

    if (!read_text_file(path, &raw_content, error)) {
        return FALSE;
    }

    wrapped_content = g_strdup_printf("[homework]\n%s", raw_content);
    key_file = g_key_file_new();
    if (!g_key_file_load_from_data(key_file, wrapped_content, -1, G_KEY_FILE_NONE, error)) {
        goto cleanup;
    }

    initialize_vectors_from_config(config, data, key_file);
    ok = TRUE;

cleanup:
    g_key_file_unref(key_file);
    g_free(wrapped_content);
    g_free(raw_content);
    return ok;
}

// Initializes the vectors x and y by asking the user for input in a style similar to the NumericalSum assignment.
static gboolean initialize_from_numericalsum_style(int argc, char **argv, Config *config, DataSet *data) {
    char output_buffer[256];

    if (argc >= 3) {
        config->vector_length = (gsize)strtoull(argv[1], NULL, 10);
        config->alpha = atof(argv[2]);
    } else {
        if (!read_gsize_value("Enter vector size N: ", &config->vector_length)) {
            return FALSE;
        }
        if (!read_double_value("Enter constant a: ", &config->alpha)) {
            return FALSE;
        }
    }

    if (!read_gsize_value("Enter chunk size: ", &config->chunk_size)) {
        return FALSE;
    }
    if (!read_string_value("Enter HDF5 output file name (.h5 is optional): ", output_buffer, sizeof(output_buffer))) {
        return FALSE;
    }

    config->output_hdf5 = normalize_hdf5_filename(output_buffer);
    if (config->vector_length == 0 || config->chunk_size == 0) {
        return FALSE;
    }

    if (!allocate_dataset(config, data, NULL)) {
        return FALSE;
    }

    return initialize_vectors_from_input(config, data);
}

// Computes the reference d vector and its sum using the formula d[i] = alpha * x[i] + y[i].
static double compute_reference(const Config *config, DataSet *data) {
    double sum = 0.0;

    for (gsize i = 0; i < config->vector_length; ++i) {
        data->d_reference[i] = config->alpha * data->x[i] + data->y[i];
        sum += data->d_reference[i];
    }

    return sum;
}

// Computes the d vector and its sum in chunks, storing partial sums for each chunk.
static double compute_chunked(const Config *config, DataSet *data) {
    double total = 0.0;

    for (gsize chunk_index = 0; chunk_index < data->number_of_chunks; ++chunk_index) {
        gsize start = chunk_index * config->chunk_size;
        gsize end = start + config->chunk_size;
        double partial = 0.0;

        if (end > config->vector_length) {
            end = config->vector_length;
        }

        for (gsize i = start; i < end; ++i) {
            data->d_chunked[i] = config->alpha * data->x[i] + data->y[i];
            partial += data->d_chunked[i];
        }

        data->partial_chunk_sum[chunk_index] = partial;
        total += partial;
    }

    return total;
}

// Compares two vectors element-wise and checks if they are within a specified tolerance.
static gboolean compare_vectors(const double *lhs, const double *rhs, gsize length, double tolerance,
                                double *max_abs_diff) {
    *max_abs_diff = 0.0;

    for (gsize i = 0; i < length; ++i) {
        double diff = fabs(lhs[i] - rhs[i]);
        if (diff > *max_abs_diff) {
            *max_abs_diff = diff;
        }
        if (diff > tolerance) {
            return FALSE;
        }
    }

    return TRUE;
}

// Writes a double dataset to an HDF5 file at the specified location and name.
#if HAVE_HDF5

// Helper function to write a double array as a dataset in an HDF5 file.
static gboolean write_double_dataset(hid_t location, const gchar *name, const double *data, hsize_t length) {
    hid_t dataspace = -1;
    hid_t dataset = -1;
    hsize_t dims[1] = {length};
    gboolean ok = FALSE;

    dataspace = H5Screate_simple(1, dims, NULL);
    if (dataspace < 0) {
        goto cleanup;
    }

    dataset = H5Dcreate2(location, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset < 0) {
        goto cleanup;
    }

    if (H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0) {
        goto cleanup;
    }

    ok = TRUE;

cleanup:
    if (dataset >= 0) {
        H5Dclose(dataset);
    }
    if (dataspace >= 0) {
        H5Sclose(dataspace);
    }
    return ok;
}

// Writes the output data to an HDF5 file, including the original vectors, the computed d vectors, and the partial chunk sums.
static gboolean write_hdf5_output(const Config *config, const DataSet *data, GError **error) {
    hid_t file = -1;
    hid_t chunks_group = -1;
    gboolean ok = FALSE;

    file = H5Fcreate(config->output_hdf5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        g_set_error(error, G_FILE_ERROR, g_file_error_from_errno(errno),
                    "failed to create HDF5 file '%s'", config->output_hdf5);
        return FALSE;
    }

    chunks_group = H5Gcreate2(file, "/chunks", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (chunks_group < 0) {
        g_set_error(error, G_FILE_ERROR, g_file_error_from_errno(errno),
                    "failed to create chunks group in '%s'", config->output_hdf5);
        goto cleanup;
    }

    if (!write_double_dataset(file, "x", data->x, config->vector_length) ||
        !write_double_dataset(file, "y", data->y, config->vector_length) ||
        !write_double_dataset(file, "d_reference", data->d_reference, config->vector_length) ||
        !write_double_dataset(file, "d_chunked", data->d_chunked, config->vector_length) ||
        !write_double_dataset(file, "partial_chunk_sum", data->partial_chunk_sum, data->number_of_chunks)) {
        g_set_error(error, G_FILE_ERROR, g_file_error_from_errno(errno),
                    "failed while writing datasets to '%s'", config->output_hdf5);
        goto cleanup;
    }

    for (gsize chunk_index = 0; chunk_index < data->number_of_chunks; ++chunk_index) {
        gchar *dataset_name = NULL;
        gsize start = chunk_index * config->chunk_size;
        gsize end = start + config->chunk_size;

        if (end > config->vector_length) {
            end = config->vector_length;
        }

        dataset_name = g_strdup_printf("chunk_%04zu", chunk_index);
        if (!write_double_dataset(chunks_group, dataset_name, &data->d_chunked[start], (hsize_t)(end - start))) {
            g_set_error(error, G_FILE_ERROR, g_file_error_from_errno(errno),
                        "failed while writing chunk datasets to '%s'", config->output_hdf5);
            g_free(dataset_name);
            goto cleanup;
        }
        g_free(dataset_name);
    }

    ok = TRUE;

cleanup:
    if (chunks_group >= 0) {
        H5Gclose(chunks_group);
    }
    if (file >= 0) {
        H5Fclose(file);
    }
    return ok;
}


#else

// If HDF5 support is not available, this function will set an error indicating that HDF5 output cannot be written.
static gboolean write_hdf5_output(const Config *config, const DataSet *data, GError **error) {
    (void)config;
    (void)data;
    g_set_error(error, G_FILE_ERROR, g_file_error_from_errno(ENOSYS),
                "HDF5 support is not available in this environment. Install libhdf5-dev and rebuild.");
    return FALSE;
}
#endif


// The main function orchestrates the loading of configuration, initialization of vectors, computation of reference and chunked sums, comparison of results, and writing of HDF5 output. 
// It also handles error reporting and resource cleanup.
int main(int argc, char **argv) {
    Config config = {0};
    DataSet data = {0};
    GError *error = NULL;
    double reference_sum;
    double chunked_sum;
    double max_abs_diff = 0.0;
    double reference_time;
    double chunked_time;
    double hdf5_time;
    double t0;
    double t1;
    gboolean vectors_match;
    gboolean hdf5_written;

    if (argc >= 2 && strcmp(argv[1], "--config") == 0) {
        const gchar *config_path = (argc >= 3) ? argv[2] : "homework6.cfg";
        if (!load_config_and_vectors(config_path, &config, &data, &error)) {
            fprintf(stderr, "Configuration error: %s\n", error->message);
            g_clear_error(&error);
            free_config(&config);
            free_dataset(&data);
            return EXIT_FAILURE;
        }
    } else {
        if (!initialize_from_numericalsum_style(argc, argv, &config, &data)) {
            fprintf(stderr, "Input error: failed to read the NumericalSum-style inputs.\n");
            free_config(&config);
            free_dataset(&data);
            return EXIT_FAILURE;
        }
    }

    t0 = wall_time_seconds();
    reference_sum = compute_reference(&config, &data);
    t1 = wall_time_seconds();
    reference_time = t1 - t0;

    t0 = wall_time_seconds();
    chunked_sum = compute_chunked(&config, &data);
    t1 = wall_time_seconds();
    chunked_time = t1 - t0;
    vectors_match = compare_vectors(data.d_reference, data.d_chunked, config.vector_length, 1e-12, &max_abs_diff);

    printf("Vector length:        %zu\n", config.vector_length);
    printf("Chunk size:           %zu\n", config.chunk_size);
    printf("Number of chunks:     %zu\n", data.number_of_chunks);
    printf("Alpha:                %.6f\n", config.alpha);
    printf("Reference sum:        %.15f\n", reference_sum);
    printf("Chunked partial sum:  %.15f\n", chunked_sum);
    printf("Max absolute diff:    %.15e\n", max_abs_diff);
    printf("d arrays match:       %s\n", vectors_match ? "YES" : "NO");
    printf("Sums match:           %s\n", fabs(reference_sum - chunked_sum) <= 1e-12 ? "YES" : "NO");
    printf("Reference time:       %.9e s\n", reference_time);
    printf("Chunked time:         %.9e s\n", chunked_time);

    t0 = wall_time_seconds();
    hdf5_written = write_hdf5_output(&config, &data, &error);
    t1 = wall_time_seconds();
    hdf5_time = t1 - t0;
    if (hdf5_written) {
        printf("HDF5 output:          written to %s\n", config.output_hdf5);
        printf("HDF5 write time:      %.9e s\n", hdf5_time);
    } else {
        fprintf(stderr, "HDF5 output:          skipped (%s)\n", error->message);
        fprintf(stderr, "HDF5 write time:      %.9e s\n", hdf5_time);
        g_clear_error(&error);
    }

    free_dataset(&data);
    free_config(&config);

    return (vectors_match && fabs(reference_sum - chunked_sum) <= 1e-12) ? EXIT_SUCCESS : EXIT_FAILURE;
}
