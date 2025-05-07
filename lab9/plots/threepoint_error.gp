# Set terminal and output
set terminal pngcairo size 800,600 enhanced font 'Arial,10'
set output 'threepoint_error.png'

# Set up the plot
set title "Threepoint method error"
set xlabel "log10(h)"
set ylabel "log10(error)"
set grid

# Plot function and data (space-separated file)
plot 'threepoint_error.csv' using 2:1 with points pt 7 ps 0.5 lc rgb "red" title "error"

# Reset output
unset output
set terminal pop
