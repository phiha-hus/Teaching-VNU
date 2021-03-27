% Shiming Duan, Aug.2, 2010
% Calculate the difference between two sides of the equation to find Qk and
% Sk
% This is the objective function used in the search for Qk
function F = LambQ(x,A,Ad,h,branch_id,reduced,notfixed)
if reduced ==0 % full elements for Q
    x1 =x;
else % Reduced elements for Q
    x1 = zeros(1,length(A)^2); 
    position = find(notfixed==1);
    for k=1:length(position)
        x1(position(k)) = x(k);
    end
    x1  = reshape(x1,length(A),[]);
end      
H = Ad*h*x1;
W = lambertw_matrix(branch_id,H);
F = W*expm(W+A*h)-Ad*h; % difference between two sides of the equation
