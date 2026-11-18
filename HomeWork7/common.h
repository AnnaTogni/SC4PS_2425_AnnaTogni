// common.h -- random number generators and helpers shared by the HomeWork7 programs
//
// As required by the exercises, the ONLY random source is the linear congruential
// generator (LCG) of the lecture notebook:
//
//     x_{n+1} = (a x_n + c) mod 2^32,   a = 1664525,  c = 1013904223
//     u_n     = x_n / 2^32  in [0, 1)
//
// The "mod 2^32" is free: uint32_t arithmetic in C wraps around modulo 2^32.
// Period: 2^32 ~ 4.3e9 (the full period, since c is odd and a - 1 is a multiple of 4).

#ifndef COMMON_H
#define COMMON_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#define OUTPUT_DIR "output"

typedef struct {
    uint32_t state;
} lcg_t;

static inline lcg_t lcg_new(uint32_t seed) { return (lcg_t){seed}; }

static inline uint32_t lcg_next(lcg_t *g)
{
    g->state = 1664525u * g->state + 1013904223u;
    return g->state;
}

// Uniform double in [0, 1): 2^32 possible values, spaced by 2^-32 ~ 2.3e-10.
static inline double lcg_uniform(lcg_t *g) { return lcg_next(g) / 4294967296.0; }

// Read an optional positive integer argument argv[i] (accepts 1e6).
// Returns def if i >= argc. Exits with a message on invalid input.
unsigned long arg_ulong(int argc, char *argv[], int i, unsigned long def, const char *name);

// Same, for a 32-bit seed (any value 0 ... 2^32 - 1).
uint32_t arg_seed(int argc, char *argv[], int i, uint32_t def);

// Create OUTPUT_DIR if needed and open OUTPUT_DIR/name for writing (exits on failure).
FILE *open_output(const char *name);

#endif
