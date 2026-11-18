// common.c -- helpers shared by the HomeWork7 programs (see common.h)

#define _POSIX_C_SOURCE 200809L // mkdir() with -std=c11

#include "common.h"

#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>

unsigned long arg_ulong(int argc, char *argv[], int i, unsigned long def, const char *name)
{
    if (i >= argc)
        return def;
    char *end;
    double v = strtod(argv[i], &end);
    if (end == argv[i] || *end != '\0' || v < 1 || v != floor(v) || v > 4e18) {
        fprintf(stderr, "Error: %s must be a positive integer, got '%s'\n", name, argv[i]);
        exit(1);
    }
    return (unsigned long)v;
}

uint32_t arg_seed(int argc, char *argv[], int i, uint32_t def)
{
    if (i >= argc)
        return def;
    char *end;
    double v = strtod(argv[i], &end);
    if (end == argv[i] || *end != '\0' || v < 0 || v != floor(v) || v > 4294967295.0) {
        fprintf(stderr, "Error: the seed must be an integer 0 ... 4294967295, got '%s'\n",
                argv[i]);
        exit(1);
    }
    return (uint32_t)v;
}

FILE *open_output(const char *name)
{
    if (mkdir(OUTPUT_DIR, 0775) != 0 && errno != EEXIST) {
        perror(OUTPUT_DIR);
        exit(1);
    }
    char path[256];
    snprintf(path, sizeof(path), "%s/%s", OUTPUT_DIR, name);
    FILE *f = fopen(path, "w");
    if (f == NULL) {
        perror(path);
        exit(1);
    }
    return f;
}
