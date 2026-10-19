#
# Tutorial Makefile for a small scientific-computing project.
#
# The purpose of this file is not only to automate compilation, but also to
# teach a few build-system ideas that are useful in research code:
#
# 1. separate source files from generated files
# 2. compile shared code once and reuse it
# 3. keep plotting reproducible with explicit targets
# 4. make external-library usage visible in one place
#

# Compiler and compiler flags.
CC := gcc
CSTD := -std=c11
WARNINGS := -Wall -Wextra -Wpedantic
OPTIMIZATION := -O2
CPPFLAGS := -Iinclude
CFLAGS := $(CSTD) $(WARNINGS) $(OPTIMIZATION)

# Query the external GNU GSL library for the include and linker flags.
# Using gsl-config is often the simplest portable approach on teaching systems.
GSL_CFLAGS := $(shell gsl-config --cflags 2>/dev/null)
GSL_LIBS := $(shell gsl-config --libs 2>/dev/null)

ifeq ($(strip $(GSL_LIBS)),)
$(error GNU GSL was not found. Install GSL and make sure gsl-config is in PATH)
endif

CPPFLAGS += $(GSL_CFLAGS)

# Keep generated artifacts separate from source files.
BUILD_DIR := build
OUTPUT_DIR := output
PLOT_DIR := plots
MPLCONFIG_DIR := .mplconfig

# COMMON_OBJ is reused by both executables.
COMMON_OBJ := $(BUILD_DIR)/fft_utils.o
PROGRAMS := $(BUILD_DIR)/sampling_demo $(BUILD_DIR)/coupled_oscillators_fft

.PHONY: all run plot plot-python clean

# Default target: build both programs.
all: $(PROGRAMS)

# Run the programs to generate numerical data in CSV form.
run: all
	./$(BUILD_DIR)/sampling_demo $(OUTPUT_DIR)
	./$(BUILD_DIR)/coupled_oscillators_fft $(OUTPUT_DIR)

# Plot with gnuplot after regenerating the data.
plot: run | $(PLOT_DIR)
	gnuplot scripts/plot_sampling.gp
	gnuplot scripts/plot_coupled_oscillators.gp

# Plot with Python after regenerating the data.
# MPLCONFIGDIR points matplotlib to a writable directory inside the project.
plot-python: run | $(PLOT_DIR) $(MPLCONFIG_DIR)
	MPLCONFIGDIR=$(MPLCONFIG_DIR) python3 scripts/plot_sampling.py
	MPLCONFIGDIR=$(MPLCONFIG_DIR) python3 scripts/plot_coupled_oscillators.py

# Directory targets.
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(OUTPUT_DIR):
	mkdir -p $(OUTPUT_DIR)

$(PLOT_DIR):
	mkdir -p $(PLOT_DIR)

$(MPLCONFIG_DIR):
	mkdir -p $(MPLCONFIG_DIR)

# Compile the shared FFT helper code once.
$(BUILD_DIR)/fft_utils.o: src/fft_utils.c include/fft_utils.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# Compile each application source file separately.
$(BUILD_DIR)/sampling_demo.o: src/sampling_demo.c include/fft_utils.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/coupled_oscillators_fft.o: src/coupled_oscillators_fft.c include/fft_utils.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# Link the executables against GSL and the C math library.
$(BUILD_DIR)/sampling_demo: $(BUILD_DIR)/sampling_demo.o $(COMMON_OBJ) | $(BUILD_DIR) $(OUTPUT_DIR)
	$(CC) $^ -o $@ $(GSL_LIBS) -lm

$(BUILD_DIR)/coupled_oscillators_fft: $(BUILD_DIR)/coupled_oscillators_fft.o $(COMMON_OBJ) | $(BUILD_DIR) $(OUTPUT_DIR)
	$(CC) $^ -o $@ $(GSL_LIBS) -lm

# Remove generated content so the project can be rebuilt from scratch.
clean:
	rm -rf $(BUILD_DIR) $(OUTPUT_DIR) $(PLOT_DIR) $(MPLCONFIG_DIR)
