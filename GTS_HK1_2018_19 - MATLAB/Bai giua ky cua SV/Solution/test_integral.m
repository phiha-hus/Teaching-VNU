function test_integral(index)
  
  switch index
    case 1
      f = @(x) x.^5;
      a = 0; b = 1;
    case 2
      f = @(x) sin(x)./x;
      a = 0; b = 1;
  end
  
  value =  Gauss_2(f,a,b)