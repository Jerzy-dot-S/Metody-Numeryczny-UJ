set term pngcairo size 1200,800
set grid
set key left top

set output 'plot_best.png'
set title 'Interpolation of f(x) = 1/(1+x^2) on [-5,5] (best N)'
set xlabel 'x'
set ylabel 'y'
plot \
    'out_best.dat' using 1:2 with lines lw 8 lc rgb '#000000' title 'f(x)', \
    'out_best.dat' using 1:3 with lines lw 2 lc rgb '#d62728' title 'Interpolacja (równoodległe)', \
    'out_best.dat' using 1:4 with lines lw 2 lc rgb '#ffffff' title 'Interpolacja (cos)', \
    'nodes_rowno.dat' using 1:2 with points pt 7 lc rgb '#d62728' title 'Węzły równoodległe', \
    'nodes_cos.dat' using 1:2 with points pt 7 lc rgb '#ffffff' title 'Węzły cosinusowe'

set output 'plot_errors.png'
set title 'Max error vs N'
set xlabel 'N (degree)'
set ylabel 'max |f - fapprox|'
set format y '%.2e'
plot \
    'errors.dat' using 1:2 with linespoints lw 2 lc rgb '#d62728' title 'Równoodległe', \
    'errors.dat' using 1:3 with linespoints lw 2 lc rgb '#2ca02c' title 'Cos'

