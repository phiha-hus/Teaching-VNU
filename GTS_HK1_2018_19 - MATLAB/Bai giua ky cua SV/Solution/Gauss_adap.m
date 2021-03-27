function [value, depth] = Gauss_adap(f,a,b,eps)
  
  if nargin == 3
    eps = 1e-7;
  end
  
  depth = 0; 
  k = 1;                    % k is the number of points between a & b
  I = Gauss_2(f,a,b);
  
  max_depth = 10;
  
  tic 

  for i = 0:max_depth-1
    depth = depth + 1 
    k = k * 2;
    h = (b-a)/k;
    diff = zeros(1,k);
    
    for j = 1:k
    temp = a + j * h - h/2;
    diff(j) = abs(Gauss_2(f,a+(j-1)*h,temp) +  Gauss_2(f,temp,a+j*h) - Gauss_2(f,a+(j-1)*h,a+j*h));
    end
  
    err = sum(diff);
    if err <= eps
      disp('Gauss adaptive method succeeded. ')
      
      Inew = 0;
      for j = 1:k
        Inew = Inew + Gauss_2(f,a+(j-1)*h,a+j*h);
      end
      value = Inew;      
      break
    else
      disp('Increase the depth by 1. ')
      i = i+1;
    end   
  end  
  
  toc