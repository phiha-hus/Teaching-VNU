reset

set term png size 1200,600
set output "poisson.png" 

set multiplot layout 1,2

set xrange [0:1]
set yrange [0:1]
set zrange [-1:1]
set cbrange [-1:1]

set pm3d; set palette

h = 1.0 / (100.0*100.0)

set title 'Original grid function'
splot "../data/gridfunction_u.dat" u 1:2:3 ti '' w l

set title 'solved grid function'
splot "../data/gridfunction_u_solve.dat" u 1:2:3 ti '' w l

unset multiplot

