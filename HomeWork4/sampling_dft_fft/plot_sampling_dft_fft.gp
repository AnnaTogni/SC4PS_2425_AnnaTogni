# plot_sampling_dft_fft.gp
#
# Reads the CSV files written by ./sampling_dft_fft in output/ and draws the same
# 2 x 3 panels as the old ROOT macro, into plots/sampling_dft_fft.png:
#
#   left column:  input signal and samples | DFT reconstruction | FFT reconstruction
#   right column: input vs observed freqs  | DFT spectrum       | FFT spectrum
#
# All panels have the same size and margins, and each column shares one x range
# (time on the left, frequency on the right), so the panels can be compared directly.

set datafile separator ","
set datafile missing "nan"
set term pngcairo size 1500,1200 enhanced font "Helvetica,10"
system "mkdir -p plots"
set output "plots/sampling_dft_fft.png"

# ---- Ranges, read from the data ----
stats "output/frequencies.csv" using 3 nooutput
f_nyq = STATS_max                                   # Nyquist frequency
stats "output/frequencies.csv" using 1 nooutput
f_max_input = STATS_max                             # largest input frequency
stats "output/continuous.csv" using 1 nooutput
t_max = STATS_max                                   # duration T
stats "output/spectrum.csv" using 3 nooutput
amp_top = (STATS_max > 0 ? 1.25 * STATS_max : 1.0)  # space above the spectral peaks
stats "output/samples.csv" using 2 nooutput
y_lim = 1.2 * (abs(STATS_min) > abs(STATS_max) ? abs(STATS_min) : abs(STATS_max))
y_lim = (y_lim > 0 ? y_lim : 1.0)

# One frequency range for the whole right column: up to the largest input
# frequency or the Nyquist frequency, whichever is larger (+5%).
f_top = 1.05 * (f_max_input > f_nyq ? f_max_input : f_nyq)

# ---- Common style ----
set grid
# "width 2" makes the legend box wider than gnuplot's own estimate of the text,
# which is too small for enhanced text and would cut it.
set key top right box opaque width 2 height 0.5 samplen 2.5 spacing 1.2
set xtics offset 0,0.2
set xlabel offset 0,0.5

# Nyquist marker, drawn on every frequency panel
set style line 10 dt 3 lw 2 lc rgb "#27ae60"

# Fixed grid: every panel gets exactly the same size and distance from its
# neighbours (margins = left,right,bottom,top of the whole figure; spacing = gaps).
set multiplot layout 3,2 title "Sampling, DFT, FFT and aliasing" \
    margins 0.06,0.98,0.05,0.93 spacing 0.08,0.09

# ======== Row 1 ========

# ---- 1. Input signal and samples ----
set title "Input signal in the time domain"
set xlabel "time [s]"
set ylabel "signal"
set xrange [0:t_max]
set yrange [-y_lim:y_lim]
plot "output/continuous.csv" using 1:2 with lines lw 1.5 lc rgb "#1f4fbf" title "continuous x(t)", \
     "output/samples.csv"    using 1:2 with points pt 7 ps 0.6 lc rgb "#c0392b" title "samples"

# ---- 2. Input frequencies vs where they are observed ----
set title "Input frequencies and where they are observed after sampling"
set xlabel "frequency [Hz]"
set ylabel "relative amplitude"
set xrange [0:f_top]
set yrange [0:1.4]
set arrow 1 from f_nyq, graph 0 to f_nyq, graph 1 nohead ls 10
set label 1 "Nyquist" at f_nyq, graph 0.5 center rotate by 90 offset -1,0 tc rgb "#27ae60"
plot "output/frequencies.csv" using 1:(1.0)  with impulses lw 3 lc rgb "#1f4fbf" title "input", \
     ""                       using 2:(0.85) with impulses lw 3 dt 2 lc rgb "#c0392b" title "observed (alias)"

# ======== Row 2 ========

# ---- 3. DFT reconstruction ----
set title "Samples reconstructed with the inverse DFT"
set xlabel "time [s]"
set ylabel "signal"
set xrange [0:t_max]
set yrange [-y_lim:y_lim]
unset arrow 1
unset label 1
plot "output/samples.csv" using 1:2 with points pt 7 ps 0.6 lc rgb "#7f8c8d" title "samples", \
     ""                   using 1:3 with lines lw 1.5 lc rgb "#8e44ad" title "IDFT(DFT(x))"

# ---- 4. DFT spectrum ----
set title "Amplitude spectrum: naive DFT"
set xlabel "frequency [Hz]"
set ylabel "amplitude"
set xrange [0:f_top]
set yrange [0:amp_top]
set arrow 1 from f_nyq, graph 0 to f_nyq, graph 1 nohead ls 10
set label 1 "Nyquist" at f_nyq, graph 0.5 center rotate by 90 offset -1,0 tc rgb "#27ae60"
plot "output/spectrum.csv" using 1:2 with impulses lw 2 lc rgb "#8e44ad" title "DFT"

# ======== Row 3 ========

# ---- 5. FFT reconstruction ----
set title "Samples reconstructed with the inverse FFT"
set xlabel "time [s]"
set ylabel "signal"
set xrange [0:t_max]
set yrange [-y_lim:y_lim]
unset arrow 1
unset label 1
plot "output/samples.csv" using 1:2 with points pt 7 ps 0.6 lc rgb "#7f8c8d" title "samples", \
     ""                   using 1:4 with lines lw 1.5 lc rgb "#e67e22" title "IFFT(FFT(x))"

# ---- 6. FFT spectrum ----
set title "Amplitude spectrum: FFT"
set xlabel "frequency [Hz]"
set ylabel "amplitude"
set xrange [0:f_top]
set yrange [0:amp_top]
set arrow 1 from f_nyq, graph 0 to f_nyq, graph 1 nohead ls 10
set label 1 "Nyquist" at f_nyq, graph 0.5 center rotate by 90 offset -1,0 tc rgb "#27ae60"
plot "output/spectrum.csv" using 1:3 with impulses lw 2 lc rgb "#e67e22" title "FFT"

unset multiplot
set output
