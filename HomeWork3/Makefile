CC     = gcc
CFLAGS = -Wall -Wextra -O3
LDLIBS = -lm

# "make" builds the program and then prints how to use it
all: matrixmult help

matrixmult: MatrixMult.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

# "make run" builds (if needed) and runs an example
run: matrixmult
	./matrixmult 2 3 1000 C.dat

# "make help" prints the instructions (the @ hides the echo command itself)
help:
	@echo ""
	@echo "Usage:   ./matrixmult a b N fileout [ones|identity]"
	@echo ""
	@echo "  a, b      value of the elements of A and B"
	@echo "  N         matrix dimension (N x N)"
	@echo "  fileout   file where the result C is saved"
	@echo "  ones      (default) every element of A is a, of B is b  ->  C_ij = N a b"
	@echo "  identity  A = a I, B = b I                              ->  C = a b I"
	@echo ""
	@echo "It computes C = A B with all 6 loop orderings, times each one,"
	@echo "checks the result and saves C to fileout."
	@echo ""
	@echo "Examples:"
	@echo "  ./matrixmult 2 3 4 C.dat              small, prints A, B and C"
	@echo "  ./matrixmult 2 3 1000 C.dat           benchmark (~20 s)"
	@echo "  ./matrixmult 2 3 1000 C.dat identity  identity matrices"
	@echo ""
	@echo "Other commands:"
	@echo "  make run     run the N = 1000 example"
	@echo "  make clean   delete the program and the output file C.dat"
	@echo ""

clean:
	rm -f matrixmult C.dat

.PHONY: all run help clean
