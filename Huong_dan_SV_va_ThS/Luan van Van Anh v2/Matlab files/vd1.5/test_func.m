h = 0.1; kp = 0.4451; K = 1.53; 
ki = 2.3046; tau = 0.0254;

A = - [ 0 1; 0 -1/tau];
Ad = - [0 0; -K*ki/tau -K*kp/tau] ; 

[alpha] = find_alpha(h,A,Ad) ; 

% K1 = J1(h) ;

%[K3] = find_K3(h, A, Ad)