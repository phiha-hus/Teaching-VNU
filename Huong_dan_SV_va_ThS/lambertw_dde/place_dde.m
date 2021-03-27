% Written by Shiming Duan, Aug. 06, 2010
% Dept. Mechanical Engineering, University of Michigan  
% Email: duansm@umich.edu                               
% For Matlab Toolbox of Time-Delay System Supplements   
% Find a feedback gain K to place the rightmost eigenvalue of the closed-loop at desired place (pole_desired)
% contr mode 1: u = -K*x(t)
% the closed-loop system is x_dot(t) = (A-B*K)*x(t) + Ad*x(t-h)
% contr mode 2: u =  -Kd*x(t-h)
% the closed-loop system is x_dot(t) = A*x(t) + (Ad-B*Kd)*x(t-h)
% contr mode 3: u = -K*x(t) - Kd*x(t-h)
% the closed-loop system is x_dot(t) = (A-B*K)*x(t) + (Ad-B*Kd)*x(t-h)
% observer mode 4: observer gain is L
% the closed-loop system is x_dot(t) = (A-L*C)*x(t) + Ad*x(t-h)
% default mode is 1.
% K = place_dde(A,Ad,B,h,branch,pole_desired,contr_mode,K0,Q_ini,thredhold)
% The arguments contr_mode, K0, Q_ini and thredhold are optional
% pole_desired must be a scalar. For complex values, e.g., you can use
% either -1+i or -1-i. They are the same here.


function K = place_dde(A,Ad,B,h,branch,pole_desired,contr_mode,K0,Q_ini,thredhold)
options_K_def = optimset('Display','final');
th2 = 0;
[B_row, B_col] = size(B);
if nargin <=5
   warning('No enough input arguments.')
elseif nargin ==6
    contr_mode = 1;% default mode 1:u(t) =  -Kx(t)
    K0 = ones(B_col,length(A));
    Q_ini = expm(-A*h); % use e^(-Ah) to initialize Q
elseif nargin ==7 %only contr mode is selected
    if contr_mode <=2
        K0 = ones(B_col,length(A)); Q_ini = expm(-A*h); % mode 1 and 2
    elseif contr_mode ==3
        K0 = [ones(B_col,length(A)) zeros(B_col,length(A))];  Q_ini = expm(-A*h); % mode 3
    elseif contr_mode ==4
        K0 = ones(length(A),B_row); Q_ini = expm(-A*h); % observer mode
    else
        warning('mode must be <= 4');
    end
elseif nargin ==8
    Q_ini = expm(-A*h);    
elseif nargin ==10
        if thredhold(1)~=0
            options_K_def.TolFun = thredhold(1);
        end
        if thredhold(2)~=0
            th2 = thredhold(2); % if thredhold is not zero, use it  to replace the default thredhold 
        end
end

if length(pole_desired)>1
    warning('Pole_desired is a scalar, please enter either a real number, or a complex number');
end


K =  fsolve(@find_K,K0,options_K_def,A,Ad,B,h,branch,pole_desired,contr_mode,Q_ini,th2) % find K



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rtn=find_K(X,A,Ad,B,h,branch,pole_desired,contr_mode,Q_initial,thredhold2)
[B_row, B_col] = size(B);
% calculate the closed-loop system
if contr_mode == 1
    K = X; Kd = zeros(B_col,length(A));A_cl = A-B*K;Ad_cl = Ad-B*Kd;
elseif contr_mode ==2
    Kd = X; K = zeros(B_col,length(A));A_cl = A-B*K;Ad_cl = Ad-B*Kd;
elseif contr_mode ==3
    K = X(:,1:length(A)); Kd = X(:,length(A)+1:end);A_cl = A-B*K;Ad_cl = Ad-B*Kd;
elseif contr_mode ==4
    A_cl = A-X*B;Ad_cl = Ad;
end

% can do the assignment for only one branch!
if length(branch)~=1
    error('Please select one and only one branch');
end

% find the solution for the closed-loop system
DDE_Sol = find_Sk(A_cl,Ad_cl,h,branch,Q_initial,thredhold2);

% showing the current result
eig_S   =  eig(DDE_Sol{1}.Sk)

% calcuate the errors of eigenvalues
if isreal(pole_desired)
    rtn = (max(real(eig_S))-pole_desired)^2; % match only the real part
else
    real_eig_diff = max(real(eig_S))-real(pole_desired);
    imag_eig_diff = max(imag(eig_S(find(real(eig_S)==max(real(eig_S))))))-abs(imag(pole_desired));
    rtn = [real_eig_diff^2;imag_eig_diff^2]; % match both real and imag parts
end    

