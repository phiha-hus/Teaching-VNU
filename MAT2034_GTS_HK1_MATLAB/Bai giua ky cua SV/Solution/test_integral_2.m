function test_integral_2(index)
  
  switch index
    case 1
      f = @(x) 1./sqrt(sin(x));
      a = 0; b = 1;
    case 2
      f = @(x) exp(1./x.^2);
      a = 0; b = 1e+6;
    case 3
      f = @(x) x .* abs(sin(1./x));
      a = 0; b = 1;
    case 4
      f = @(x) x.^5;
      a = 0; b = 1;
  end
  
 %value =  Gauss_2(f,a,b)
 
 eps = 1e-10;
 [value,depth] =  Gauss_adap(f,a,b,eps)