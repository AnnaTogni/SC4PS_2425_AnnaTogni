CC     = gcc
CFLAGS = -Wall -Wextra -O2
LDLIBS = -lm

# "make" builds the program and then prints how to use it
all: numericalsum help

numericalsum: NumericalSum.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

# "make help" prints the instructions (the @ hides the echo command itself)
help:
	@echo ""
	@echo "Usage:   ./numericalsum N a x y [expected]"
	@echo ""
	@echo "  N         vector dimension (e.g. 10, 1e6, 1e8)"
	@echo "  a         scalar"
	@echo "  x, y      value of every element of x and y"
	@echo "  expected  (optional) value every d[i] should be equal to"
	@echo ""
	@echo "Examples:"
	@echo "  ./numericalsum 10  3 0.1 7.1 7.4"
	@echo "  ./numericalsum 1e8 3 0.1 7.1 7.4"
	@echo ""
	@echo "Other commands:  make clean  (delete the program)"
	@echo ""

clean:
	rm -f numericalsum

.PHONY: all help clean
