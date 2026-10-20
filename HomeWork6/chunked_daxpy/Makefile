CC     = gcc
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -O2
# GLib (config-file parser) and HDF5: pkg-config gives the right paths
PKGS   = glib-2.0 hdf5
CPPFLAGS = $(shell pkg-config --cflags $(PKGS))
LDLIBS   = $(shell pkg-config --libs $(PKGS)) -lm

# "make" builds the program and then prints how to use it
all: chunked_daxpy help

chunked_daxpy: chunked_daxpy.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(LDLIBS)

# "make run" uses the example input file daxpy.cfg
run: chunked_daxpy
	./chunked_daxpy daxpy.cfg

# "make test" runs the tests in tests/run_tests.sh
test: chunked_daxpy
	bash tests/run_tests.sh

# "make help" prints the instructions (the @ hides the echo command itself)
help:
	@echo ""
	@echo "Usage:   ./chunked_daxpy [config_file]      (default: daxpy.cfg)"
	@echo ""
	@echo "The config file has one 'Variable = Value' per line:"
	@echo "  n                vector length"
	@echo "  chunk_size       elements per chunk        (give this OR number_of_chunks)"
	@echo "  number_of_chunks how many chunks"
	@echo "  a, x, y          d[i] = a*x[i] + y[i], x[i] = x + i*x_step, y[i] = y + i*y_step"
	@echo "  x_step, y_step   (optional, default 0)"
	@echo "  output           HDF5 file name           (optional, default daxpy_chunks.h5)"
	@echo ""
	@echo "It computes d with one loop and chunk by chunk, checks they are identical,"
	@echo "compares the sum of the partial chunk sums with the single-loop sum, and saves"
	@echo "everything in the HDF5 file (look inside with: h5ls -r file.h5, h5dump file.h5)."
	@echo ""
	@echo "Other commands:"
	@echo "  make run     run with daxpy.cfg"
	@echo "  make test    run the tests (edge cases and wrong input files)"
	@echo "  make clean   delete the program and the .h5 files"
	@echo ""

clean:
	rm -f chunked_daxpy *.h5
	rm -rf tests/tmp

.PHONY: all run test help clean
