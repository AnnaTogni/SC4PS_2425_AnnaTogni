CC     = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2
LDLIBS = -lgsl -lgslcblas -lm

# "make" builds the program and then prints how to use it
all: floating_point_sum help

floating_point_sum: floating_point_sum.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

# "make run" sums the lecture's vector {1, 1e16, -1e16, -0.5}
run: floating_point_sum
	./floating_point_sum

# "make help" prints the instructions (the @ hides the echo command itself)
help:
	@echo ""
	@echo "Usage:   ./floating_point_sum [x1 x2 ...]"
	@echo ""
	@echo "  without arguments: sums the lecture's vector  1 1e16 -1e16 -0.5  (exact sum 0.5)"
	@echo "  with arguments:    sums the numbers you give"
	@echo ""
	@echo "Compares 8 ways of summing: for loop, GSL, Kahan, Neumaier, pairwise,"
	@echo "sorted by |x|, long double, and the exact (correctly rounded) sum."
	@echo "For short vectors it also prints Kahan and Neumaier step by step."
	@echo ""
	@echo "Examples:"
	@echo "  ./floating_point_sum                    lecture example"
	@echo "  ./floating_point_sum 1e16 1 -1e16       the 1 is lost by the plain loop"
	@echo "  ./floating_point_sum 0.1 0.2 0.3 -0.6   not 0: 0.1, 0.2 ... are not exact in binary"
	@echo ""
	@echo "Other commands:"
	@echo "  make run     run the lecture example"
	@echo "  make clean   delete the program"
	@echo ""

clean:
	rm -f floating_point_sum

.PHONY: all run help clean
