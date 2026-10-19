# plot_coupled_oscillators.gp
#
# Tutorial plotting script for the coupled-oscillator example.
#
# It generates:
# 1. the displacement as a function of time
# 2. the FFT amplitude spectrum of x1(t)
#
# Compared with the Python version, this script is shorter and less flexible,
# but it is useful as a lightweight plotting option on scientific Linux
# systems.

set datafile separator ","
set grid
set key outside
set term pngcairo size 1400,700 enhanced font "Helvetica,11"

system "mkdir -p plots"

set output "plots/coupled_oscillators_time.png"
set multiplot layout 2,1 rowsfirst title "Coupled oscillators: displacement versus time"

# First panel: both masses over the full simulation time.
set xlabel "Time (s)"
set ylabel "Displacement (m)"
set title "Mass 1 and mass 2"
plot "output/coupled_oscillators_time.csv" using 1:2 with lines lw 2 title "x1(t)", \
     "output/coupled_oscillators_time.csv" using 1:3 with lines lw 2 title "x2(t)"

# Second panel: a short-time window is easier to read at a glance.
set title "Mass 1 over an early-time window"
plot [0:20] "output/coupled_oscillators_time.csv" using 1:2 with lines lw 2 title "x1(t)"

unset multiplot

set output "plots/coupled_oscillators_spectrum.png"
set xlabel "Frequency (Hz)"
set ylabel "Amplitude"
# The Python version also overlays the theoretical mode frequencies.
set title "Coupled oscillators: FFT amplitude spectrum of x1(t)"
plot [0:2] "output/coupled_oscillators_spectrum.csv" using 1:2 with lines lw 2 title "FFT amplitude"

set output
