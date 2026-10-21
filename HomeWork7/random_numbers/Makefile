CC     = gcc
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -O2
LDLIBS = -lgsl -lgslcblas -lm

PROGRAMS = coin_tosses monte_carlo_pi square_uniform exponential empirical_cdf

# "make" builds the five programs and then prints how to use them
all: $(PROGRAMS) help

# every program = its own .c file + the shared helpers in common.c
$(PROGRAMS): %: %.c common.c common.h
	$(CC) $(CFLAGS) -o $@ $< common.c $(LDLIBS)

exponential empirical_cdf: exponential.h

# "make run" runs the five exercises with the default parameters
run: $(PROGRAMS)
	./coin_tosses
	./monte_carlo_pi
	./square_uniform
	./exponential
	./empirical_cdf

# "make plot" runs everything and draws the figures in plots/
plot: run
	MPLCONFIGDIR=/tmp/mplconfig python3 plot_random.py

# "make help" prints the instructions (the @ hides the echo command itself)
help:
	@echo ""
	@echo "HomeWork7: random numbers. Every program accepts optional arguments;"
	@echo "the last one is always the seed of the generator (default 12345)."
	@echo ""
	@echo "  ./coin_tosses    [N] [seed]          1. N coin tosses, running fraction of heads"
	@echo "  ./monte_carlo_pi [N_max] [R] [seed]  2. pi from N = 10 ... N_max points, R runs each"
	@echo "  ./square_uniform [n] [bins] [seed]   3. Y = U^2: histogram vs 1/(2 sqrt(y))"
	@echo "  ./exponential    [n] [bins] [seed]   4. Y = -ln(1-U)/1.5: histogram vs the pdf"
	@echo "  ./empirical_cdf  [n] [seed]          5. empirical CDF of the sample of 4"
	@echo ""
	@echo "Defaults: N = 1e5, N_max = 1e7, R = 20, n = 1e6. Results go to output/*.csv."
	@echo ""
	@echo "Other commands:"
	@echo "  make run     run the five exercises"
	@echo "  make plot    run them and draw plots/*.png (needs Python with matplotlib)"
	@echo "  make clean   delete the programs and output/ (plots/ is kept)"
	@echo ""

clean:
	rm -f $(PROGRAMS)
	rm -rf output

.PHONY: all run plot help clean
