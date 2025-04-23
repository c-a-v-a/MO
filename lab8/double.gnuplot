set terminal pngcairo size 800,600 enhanced
set output 'img/double.png'

set title "Double"
set xlabel "log10(h)"
set ylabel "log10(|err|)"

plot 'double.dat' using 1:2 with linespoints title 'x=0 forward', \
		'' using 1:3 with linespoints title 'x=0 threepoint', \
		'' using 1:4 with linespoints title 'x=pi/4 forward', \
		'' using 1:5 with linespoints title 'x=pi/4 backward', \
		'' using 1:6 with linespoints title 'x=pi/4 central', \
		'' using 1:7 with linespoints title 'x=pi/2 backward', \
		'' using 1:8 with linespoints title 'x=pi/2 threepoint'
