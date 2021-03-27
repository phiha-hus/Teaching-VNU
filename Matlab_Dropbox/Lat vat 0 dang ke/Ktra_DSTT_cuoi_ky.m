% De kiem tra dai so tuyen tinh cuoi ky

syms x

 A = [1 x 3 2;0 4 -2 -3;0 0 2*x-1 x;0 0  0 3*x]
 
  B = [1 1 1 1;0 1 0 -1;0 0 -1 0;0 0 0 -1]'
  
   C=B*A
   
   det(C)
   
   expand(ans)
   
   % solve the equation det(C)=10x-4