CC = gcc
CFLAGS = -O2 -Wall -Wextra -Wpedantic -std=c11
LDLIBS = -lgsl -lgslcblas -lquadmath -lm

TARGET = legendre_stability
SRC = legendre_stability.c
PYTHON_PLOT_SCRIPT = plot_legendre.py

.PHONY: all run plot python_plot python pplot clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LDLIBS)

run: $(TARGET)
	./$(TARGET)

plot: python_plot

python_plot: run $(PYTHON_PLOT_SCRIPT)
	mkdir -p plots_python
	MPLCONFIGDIR=/tmp/mplconfig python3 $(PYTHON_PLOT_SCRIPT)

python: python_plot

pplot: python_plot

clean:
	rm -f $(TARGET) legendre_x_*.csv
	rm -rf plots_python
