reset

set term png size 1000,1000
set output "euler_rk_classic.png"
set xrange [0:20]

set multiplot layout 2,2

set xlabel 'time / s'
set ylabel 'position x(t)'
set yrange [-6:6]
set grid

plot "../data/ode_euler_delta0.4000.dat" u 1:2 w l ti 'euler, delta = 0.4000', \
     "../data/ode_euler_delta0.1000.dat" u 1:2 w l ti 'euler, delta = 0.1000', \
     "../data/ode_euler_delta0.0250.dat" u 1:2 w l ti 'euler, delta = 0.0250', \
     "../data/ode_euler_delta0.0063.dat" u 1:2 w l ti 'euler, delta = 0.0064', \
     "../data/ode_euler_delta0.0016.dat" u 1:2 w l ti 'euler, delta = 0.0016'  
     
set xlabel 'time / s'
set ylabel 'position x(t)'
set yrange [-6:6]
set grid

plot "../data/ode_rk_classic_delta0.4000.dat" u 1:2 w l ti 'rk classic, delta = 0.4000', \
     "../data/ode_rk_classic_delta0.1000.dat" u 1:2 w l ti 'rk classic, delta = 0.1000', \
     "../data/ode_rk_classic_delta0.0250.dat" u 1:2 w l ti 'rk classic, delta = 0.0250', \
     "../data/ode_rk_classic_delta0.0063.dat" u 1:2 w l ti 'rk classic, delta = 0.0064', \
     "../data/ode_rk_classic_delta0.0016.dat" u 1:2 w l ti 'rk classic, delta = 0.0016'  

set xlabel 'time / s'
set ylabel 'velocity v(t)'
set yrange [-6:6]
set grid 

plot "../data/ode_euler_delta0.4000.dat" u 1:3 w l ti 'euler, delta = 0.4000', \
     "../data/ode_euler_delta0.1000.dat" u 1:3 w l ti 'euler, delta = 0.1000', \
     "../data/ode_euler_delta0.0250.dat" u 1:3 w l ti 'euler, delta = 0.0250', \
     "../data/ode_euler_delta0.0063.dat" u 1:3 w l ti 'euler, delta = 0.0064', \
     "../data/ode_euler_delta0.0016.dat" u 1:3 w l ti 'euler, delta = 0.0016'
     
set xlabel 'time / s'
set ylabel 'velocity v(t)'
set yrange [-6:6]
set grid 

plot "../data/ode_rk_classic_delta0.4000.dat" u 1:3 w l ti 'rk classic, delta = 0.4000', \
     "../data/ode_rk_classic_delta0.1000.dat" u 1:3 w l ti 'rk classic, delta = 0.1000', \
     "../data/ode_rk_classic_delta0.0250.dat" u 1:3 w l ti 'rk classic, delta = 0.0250', \
     "../data/ode_rk_classic_delta0.0063.dat" u 1:3 w l ti 'rk classic, delta = 0.0064', \
     "../data/ode_rk_classic_delta0.0016.dat" u 1:3 w l ti 'rk classic, delta = 0.0016'
   

unset multiplot
