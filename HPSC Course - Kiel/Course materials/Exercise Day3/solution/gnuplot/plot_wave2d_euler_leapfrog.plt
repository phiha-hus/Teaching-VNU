reset

set xrange [0:1]
set yrange [0:1]
set zrange [-1:1]
set cbrange [-1:1]

set term gif animate size 1200,600
set output "wave2d.gif" 

set pm3d; set palette

n = 50
h = 1.0 / (n*n)

do for [ii=0:10000:100] {
  t = ii * h
  set multiplot layout 1,2
  set title "Wave equation 2D, euler"
  splot sprintf("../data/wave2d_euler_%05d.dat",ii) u 1:2:3 ti sprintf("t = %05.2f s", t) w l
  set title "Wave equation 2D, leapfrog"
  splot sprintf("../data/wave2d_leapfrog_%05d.dat",ii) u 1:2:3 ti sprintf("t = %05.2f s", t) w l
  unset multiplot
}


