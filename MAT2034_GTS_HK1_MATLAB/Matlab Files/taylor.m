function a = taylor(n,x0)
  
 a = input('Input the coefficients of the input polynomial in the vector form [a0 a1 a2 ... an]: ')
  
  if length(a) ~= n+1 
    error('do dai cua vecto he so khong chinh xac')
  end
  
  for k = 1:n
    for j = n : -1 : k
      a(j) = a(j) + x0 * a(j+1);
     end     
  end
