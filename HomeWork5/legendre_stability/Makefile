CC     = gcc
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -O2
LDLIBS = -lgsl -lgslcblas -lquadmath -lm

# "make" builds the program and then prints how to use it
all: legendre_stability help

legendre_stability: legendre_stability.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

# "make run" runs the default study (CSV files in output/)
run: legendre_stability
	./legendre_stability

# "make plot" runs the default study and draws the figures in plots/
plot: run
	MPLCONFIGDIR=/tmp/mplconfig python3 plot_legendre.py

# "make help" prints the instructions (the @ hides the echo command itself)
help:
	@echo ""
	@echo "Usage:   ./legendre_stability [l_max] [L] [x1 x2 ...]"
	@echo ""
	@echo "  l_max   largest degree l                              (default 50)"
	@echo "  L       start of Miller's backward recurrence, L > l_max (default 80)"
	@echo "  x_i     |x| < 1 studies P_l(x), x > 1 studies Q_l(x)"
	@echo "          (default 0.1 0.5 0.9 0.99 1.1 1.5 3)"
	@echo ""
	@echo "For each x it compares the forward and the backward recurrence (Miller with P_0 only,"
	@echo "and, for |x| < 1, with P_0 and P_1 fixed as in the professor's correction) in double"
	@echo "precision with a quad-precision reference, prints the largest errors and writes"
	@echo "output/legendre_<P|Q>_x<value>.csv. Then: python3 plot_legendre.py"
	@echo ""
	@echo "Examples:"
	@echo "  ./legendre_stability                    the homework study (+ Q_l for x > 1)"
	@echo "  ./legendre_stability 50 200 1.1         Miller needs a larger L when x -> 1"
	@echo "  ./legendre_stability 1000 1500 0.99     forward still stable for large l"
	@echo ""
	@echo "Other commands:"
	@echo "  make run     run the default study"
	@echo "  make plot    run it and draw plots/*.png"
	@echo "  make clean   delete the program, output/ and plots/"
	@echo ""

clean:
	rm -rf legendre_stability output plots

.PHONY: all run plot help clean
