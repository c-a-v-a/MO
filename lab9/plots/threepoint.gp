# Set terminal and output
set terminal pngcairo size 800,600 enhanced font 'Arial,10'
set output 'threepoint.png'

# Define constants
sqrt5 = sqrt(5)
coth_sqrt5 = cosh(sqrt5)/sinh(sqrt5)

# Define the function U(x)
U(x) = -((9 \
 - 95*exp((-1 - sqrt5)*(-1 + x)) \
 + 55*exp((-1 + sqrt5)*x) \
 + 95*exp(1 + sqrt5 + (-1 + sqrt5)*x) \
 - 55*exp(2*sqrt5 - (1 + sqrt5)*x) \
 + 2*x*(6 + x*(3 + 2*x)) \
 - exp(2*sqrt5)*(9 + 2*x*(6 + x*(3 + 2*x))) \
)*(-1 + coth_sqrt5))/64

# Set up the plot
set title "U(x) for threepoint"
set xlabel "x"
set ylabel "U(x)"
set grid

# Plot function and data (space-separated file)
plot [0:1] U(x) title "U(x)" lw 2 lc rgb "blue", \
     'threepoint.csv' using 1:2 with points pt 7 ps 0.5 lc rgb "red" title "Approximation"

# Reset output
unset output
set terminal pop
