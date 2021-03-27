% Written by Shiming Duan, Aug. 08, 2010
% Dept. Mechanical Engineering, University of Michigan  
% Email: duansm@umich.edu                               
% For Matlab Toolbox of Time-Delay System Supplements   
% Calculate the stability radius for a perturbed time delay system
% x_dot(t) = (A+E*Delta1*F1)x(t) + (Ad+E*Delta2*F2)x(t-h)


function sr = stabilityradiu_dde(A,Ad,E,F1,F2,h)
W0=[F1;F2]*inv(-A-Ad)*E;
om_star=max(svd(A))+max(svd(Ad))+max(svd(E))*min(svd(W0))*max(svd([F1' F2']));
stepsize = om_star/100;
k=0;
options = optimset('Display','off');
for om=0:stepsize:om_star
    k=k+1;
    [x,fval] = fminbnd(@find_r,0,1,options,A,Ad,E,F1,F2,h,om); % find the maximum of largest singular value of RM over [0,1]
    Sn(k)=fval;
end
sr=1/max(Sn);

%% rstab 
function f = find_r(x,A,Ad,E,F1,F2,h,om)
s=complex(0,om);
omega=[F1;F2*exp(-h*s)]*inv(s*eye(length(A))-A-Ad*exp(-s*h))*E;
RM=[real(omega) -x*imag(omega);1/x*imag(omega) real(omega)];
S=svd(RM,0);
f =max(S);
                