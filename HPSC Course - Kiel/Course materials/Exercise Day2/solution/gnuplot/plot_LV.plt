reset

set term png size 2000,1000
set output "lotka_volterra.png"

set multiplot layout 2,3

set xlabel 'Number of prey'
set ylabel 'Number of predators'
set xrange[0:4]
set yrange[0:6]
set grid

plot "../data/LV_euler_delta0.4000.dat" u 2:3 w l ti 'Euler, delta = 0.4000', \
     "../data/LV_runge_delta0.4000.dat" u 2:3 w l ti 'Runge, delta = 0.4000', \
     "../data/LV_rk4_delta0.4000.dat" u 2:3 w l ti 'RK 4, delta = 0.4000'
     
plot "../data/LV_euler_delta0.1000.dat" u 2:3 w l ti 'Euler, delta = 0.1000', \
     "../data/LV_runge_delta0.1000.dat" u 2:3 w l ti 'Runge, delta = 0.1000', \
     "../data/LV_rk4_delta0.1000.dat" u 2:3 w l ti 'RK 4, delta = 0.1000'
     
plot "../data/LV_euler_delta0.0250.dat" u 2:3 w l ti 'Euler, delta = 0.0250', \
     "../data/LV_runge_delta0.0250.dat" u 2:3 w l ti 'Runge, delta = 0.0250', \
     "../data/LV_rk4_delta0.0250.dat" u 2:3 w l ti 'RK 4, delta = 0.0250'
     
plot "../data/LV_euler_delta0.0063.dat" u 2:3 w l ti 'Euler, delta = 0.0063', \
     "../data/LV_runge_delta0.0063.dat" u 2:3 w l ti 'Runge, delta = 0.0063', \
     "../data/LV_rk4_delta0.0063.dat" u 2:3 w l ti 'RK 4, delta = 0.0063'
     
plot "../data/LV_euler_delta0.0016.dat" u 2:3 w l ti 'Euler, delta = 0.0016', \
     "../data/LV_runge_delta0.0016.dat" u 2:3 w l ti 'Runge, delta = 0.0016', \
     "../data/LV_rk4_delta0.0016.dat" u 2:3 w l ti 'RK 4, delta = 0.0016'
 
 
unset multiplot
