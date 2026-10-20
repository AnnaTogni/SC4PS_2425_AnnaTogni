CC     = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2
LDLIBS = -lm

# "make" builds the program and then prints how to use it
all: sampling_dft_fft help

sampling_dft_fft: sampling_dft_fft.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

# "make run" runs the classic aliasing example (50 Hz + 120 Hz sampled at 128 Hz)
run: sampling_dft_fft
	./sampling_dft_fft 128 1 50 120

# "make plot" runs the example and draws the figure with gnuplot
plot: run
	gnuplot plot_sampling_dft_fft.gp
	@echo "Figure written to plots/sampling_dft_fft.png"

# "make help" prints the instructions (the @ hides the echo command itself)
help:
	@echo ""
	@echo "Usage:   ./sampling_dft_fft fs T f1 [f2 ...]"
	@echo ""
	@echo "  fs    sampling frequency [Hz]"
	@echo "  T     duration of the sampled signal [s]   (N = fs * T samples, any N)"
	@echo "  f_i   frequencies [Hz] of the signal  x(t) = sum_i sin(2 pi f_i t)"
	@echo ""
	@echo "It predicts the aliasing, computes the spectrum with a naive DFT and an FFT,"
	@echo "times them, checks that they agree and writes CSV files in output/."
	@echo ""
	@echo "Examples:"
	@echo "  ./sampling_dft_fft 128 1 50 120      120 Hz above Nyquist -> alias at 8 Hz"
	@echo "  ./sampling_dft_fft 512 1 50 120      well sampled"
	@echo "  ./sampling_dft_fft 100 1.5 20 37     N = 150, not a power of two"
	@echo "  ./sampling_dft_fft 128 1 64          sine exactly at Nyquist: invisible"
	@echo ""
	@echo "Then draw the figure of the last run:  gnuplot plot_sampling_dft_fft.gp"
	@echo ""
	@echo "Other commands:"
	@echo "  make run     run the 50 Hz + 120 Hz example"
	@echo "  make plot    run the example and draw plots/sampling_dft_fft.png"
	@echo "  make clean   delete the program, output/ and plots/"
	@echo ""

clean:
	rm -rf sampling_dft_fft output plots

.PHONY: all run plot help clean
