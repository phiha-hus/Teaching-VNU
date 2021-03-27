% The wrong procedure, it causes missing some algebraic constraints
% Phi Monday 8, Dec 2014

%E = [1 0;0 0;0 1]
%A = [1 1;0 0;1 0]

%E = [1 0;0 0;0 0]
%A = [1 1;0 0;1 0]

%E = [0 0 1;1 0 0;0 0 0]
%A = [0 0 0;0 0 0;1 0 0]
%B = [0 0 0;0 1 0;0 0 0]


E = [0 1;0 0]
A = [1 0;0 1]
B = [0 0;0 1]

%P = rand(3,3)
%Q = rand(3,3) 
%P = rand(2,2)
%Q =rand(2,2)
%E = P * E * Q
%A = P * A * Q

n = size(E,2);

M = [-A E zeros(n);zeros(n) -A E]
P = [B zeros(n,2*n); zeros(n) B zeros(n)]
% M = [-A E zeros(n,4*n);zeros(n) -A E zeros(n,3*n); B zeros(n,2*n) -A E zeros(n); zeros(n) B zeros(n,2*n) -A E]
% P = [B zeros(n,5*n); zeros(n) B zeros(n,4*n); zeros(2*n,6*n)]

tol = 1e-7

M2 = transpose(M(:,2*n+1 : end))

U1 = null(M2,tol)
U2 = orth(U1' * M(:,1:2*n),tol)
U = U1 * U2

disp('Check non-advancedness:')
if U' * P(:,n+1 : end) == 0
	disp('Non-advanced')
	else disp('May be advanced')
end

Mt = U' * M(:,(n+1):2*n)
Nt = U' * M(:,1:n)
Pt = U' * P(:,1:n)

Z2 = null(Mt',tol)
Y2 = orth(Z2' * Nt,tol)

if isnan(Y2) == false
	disp('The set of algebraic constraints is') 
                     A2 = Y2' * Z2' * Nt
                     B2 = Y2' * Z2' * Pt
                     
                     T2 = null(Z2' * Nt,tol)
                     Z1 = orth(Mt * T2,tol)
                     
                     if min(size(Z1)) == 0
      		disp('No differential equations')
      		else
      		disp('The set of differential equations is')
      		E1 = Z1' * Mt
      		A1 = Z1' * Nt
      		B1 = Z1' * Pt
                     end
      else 
      	disp('No algebraic constraints')
      	Z1 = orth(Mt,tol)
      	E1 = Z1' * Mt
      	if isnan(E1) == true
      		disp('No differential equations')
      		else
      		disp('The set of differential equations is')
      		E1
      		A1 = Z1' * Nt
      		B1 = Z1' * Pt      		
                     end
end

