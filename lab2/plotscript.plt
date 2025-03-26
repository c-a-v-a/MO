set terminal pngcairo enhanced size 800,600

set output 'plot.png'

set title "Wykres"
set xlabel "Log10(x)"
set ylabel "Log10(e)"

set grid

plot "out.txt" using 1:2 with lines title "Data"
