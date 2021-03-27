function rtn = equations(x)

global h k1 k2 A Ad

X     = [x(1) x(2) x(3) x(4)];
Q     =  [x(1) x(2);x(3) x(4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Di tim ham W_k(A_d \tau Q_k) trong phuong trinh (10))
D =  Ad*h*Q;
[v,d] =  eig(D);
w11 = lambertw(k1,d(1,1)); 
W     =  v * [w11 w11/(d(1,1) * (1+w11));0 lambertw(k2,d(2,2))] * inv(v); % for matrix Lambert W function (Eq. (2.16)) in book
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


EX    =  expm(W+A*h);
Left  =  W * EX;                                           % the left side of Eq. (2.19) in book
Right =  Ad*h;                                             % the right side of Eq. (2.19) in book 

rtn   =  [Left(1,1)-Right(1,1)
          Left(1,2)-Right(1,2)    
          Left(2,1)-Right(2,1)
          Left(2,2)-Right(2,2)];                              % Eq. (2.19) in book 
