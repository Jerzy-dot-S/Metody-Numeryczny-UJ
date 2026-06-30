set term pngcairo size 1200,800
set grid
set key left top

set output 'N07_solution.png'
set title 'Rozwiązanie układu (zadanie N07)'
set xlabel 'x'
set ylabel 'u(x)'
plot \
    'N07_solution.dat' using 1:2 with lines lw 3 lc rgb '#d62728' title 'u(x)'
