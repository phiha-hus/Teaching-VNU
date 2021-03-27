function [value] = Gauss_2(f,a,b)
  
  s = sqrt(3);
  t = [-1/s  1/s];
  
  x = (b-a)/2 * t + (b+a)/2;
  
  % vi dx = (b-a)/2 * dt
  
  value = (b-a)/2 * sum(f(x));
