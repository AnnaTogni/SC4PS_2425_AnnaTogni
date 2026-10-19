# plot_sampling.gp
#
# Tutorial plotting script for the sampling demonstration.
#
# This script reads the CSV files produced by src/sampling_demo.c and generates
# two figures:
# 1. the sampled signals in the time domain
# 2. their one-sided FFT amplitude spectra
#
# The script is intentionally explicit rather than clever. The goal is that a
# student who has never used gnuplot can still read it and understand how the
# plot is assembled.

set datafile separator ","
set grid
set key outside
set xlabel "Time (s)"
set ylabel "Signal"
set term pngcairo size 1400,900 enhanced font "Helvetica,11"

system "mkdir -p plots"

set output "plots/sampling_signals.png"
set multiplot layout 3,1 rowsfirst title "Sampling demo: time-domain signals"

# Good sampling: both true frequencies are below Nyquist.
set title "Good sampling: 50 Hz and 120 Hz resolved correctly"
plot "output/good_sampling_signal.csv" using 1:2 every 2 with lines lw 2 title "signal(t)"

# Undersampling: the 120 Hz component is aliased.
set title "Undersampled: 120 Hz aliases to about 8 Hz"
plot "output/undersampled_signal.csv" using 1:2 with lines lw 2 title "signal(t)"

# Short total duration: the nearby frequencies are not well separated.
set title "Short record: nearby frequencies cannot be separated well"
plot "output/short_record_signal.csv" using 1:2 with lines lw 2 title "signal(t)"

unset multiplot

set xlabel "Frequency (Hz)"
set ylabel "Amplitude"
set output "plots/sampling_spectra.png"
set multiplot layout 3,1 rowsfirst title "Sampling demo: FFT amplitude spectra"

# Plot the one-sided spectra with impulses to emphasize discrete FFT bins.
set title "Good sampling spectrum"
plot [0:200] "output/good_sampling_spectrum.csv" using 1:2 with impulses lw 2 title "FFT amplitude"

set title "Undersampled spectrum"
plot [0:80] "output/undersampled_spectrum.csv" using 1:2 with impulses lw 2 title "FFT amplitude"

set title "Short-record spectrum"
plot [0:120] "output/short_record_spectrum.csv" using 1:2 with impulses lw 2 title "FFT amplitude"

unset multiplot
set output
