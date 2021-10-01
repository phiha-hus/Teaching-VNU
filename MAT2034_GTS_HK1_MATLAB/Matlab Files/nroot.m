function [x,n] = nroot(m,a,x)
  
  eps = 1e-6;
  nmax = 1e+3;
  n = 0;
  
  for i=1:nmax
    xnew = x - (x.^m-a)/(m*(x.^(m-1)));
    diff = abs(xnew - x) ;
    
    if diff <= eps
      break
    else
      x = xnew ;
      n = n + 1 ;   
    end
  
  end
  