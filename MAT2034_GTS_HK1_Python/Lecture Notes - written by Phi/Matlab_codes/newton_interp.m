function value = newton_interp(x,y,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton Interpolation Using Divided Difference Scheme
% Form: value = newton_interp(x,y,w)
% x: interpolation nodes 
% y: value of the function at nodes in vector x
% w : point where we need to compute f(w) using Newton interpolation scheme 
% value: Output interpolated value of f(w), which is P(w)
% Supporting function: div_diff
% Copyright: Phi Ha, October 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = div_diff(x,y)
n = length(v);

value = v(n);

for i = n-1 : -1 : 1
  value = v(i) + (w-x(i)) * value; 
end   

  