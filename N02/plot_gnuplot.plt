set datafile separator ','

# Timing plot (log-log)
set terminal pngcairo size 1000,700 enhanced font 'Verdana,11'
set output 'times_N02.png'
set title 'Execution time vs N for N02 (log-log)'
set xlabel 'N'
set ylabel 'Time (s)'
set grid
set logscale x 10
set logscale y 10
set key left top
set datafile missing '-1'
# Skip non-positive values on log scale via conditional (1/0 -> undefined)
plot 'times_N02.csv' using 1:(($2>0)?$2:1/0) with linespoints lt 1 lw 2 pt 7 title 'Dense', \
     '' using 1:(($3>0)?$3:1/0) with linespoints lt 2 lw 2 pt 7 title 'Woodbury', \
     '' using 1:(($4>0)?$4:1/0) with linespoints lt 3 lw 2 pt 7 title 'SOR'

# Solution plot for N = 1000
set terminal pngcairo size 1000,700 enhanced font 'Verdana,11'
set output 'solution_N02_1000.png'
set title 'Solution for N = 1000 (x_n = (n-1)h)'
set xlabel 'x'
set ylabel 'u(x)'
set grid
unset logscale
set key bottom right
plot 'solution_N1000.csv' using 1:2 with lines lt 1 lw 1 title 'dense', \
     'solution_N1000.csv' using 1:3 with lines lt 2 lw 1 title 'woodbury', \
     'solution_N1000.csv' using 1:4 with lines lt 3 lw 1 title 'sor'

# End