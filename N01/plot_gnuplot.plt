set datafile separator ','

# Timing plot (log-log)
set terminal pngcairo size 1000,700 enhanced font 'Verdana,11'
set output 'times_N01.png'
set title 'Execution time vs N for N01 (log-log)'
set xlabel 'N'
set ylabel 'Time (s)'
set grid
set logscale x 10
set logscale y 10
set key left top
plot 'times_N01.csv' using 1:4 with linespoints lt 1 lw 2 pt 7 title 'Gauss-Seidel (SOR)', \
     'times_N01.csv' using 1:3 with linespoints lt 2 lw 2 pt 7 title 'Thomas', \
     'times_N01.csv' using 1:2 with linespoints lt 3 lw 2 pt 7 title 'Dense'

# Solution plot for N = 1000
set terminal pngcairo size 1000,700 enhanced font 'Verdana,11'
set output 'solution_N1000_points.png'
set title 'Solution points for N = 1000 (x_n = (n-1)h, u_n = u(x))'
set xlabel 'x'
set ylabel 'u(x)'
set grid
unset logscale
set key bottom right
plot 'solution_N1000.csv' using 1:2 with linespoints lt 4 lw 1 pt 6 ps 0.6 title 'u_n (discrete)'

# End