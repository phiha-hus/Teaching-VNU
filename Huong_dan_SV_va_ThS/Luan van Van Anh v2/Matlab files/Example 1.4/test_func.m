global	h	A	Ad

kp = 0.4451; K = 1.53; 
ki = 2.3046; tau = 0.0254;

%h = 1; A = [1 3;-2 5]; Ad = [-1.66 0.697;-0.93 0.33]; 
h = 5; A = [0 1;-5 -1]; Ad = [0 0;-3 -0.6];

%rank(Ad)

[alpha] = find_alpha(A,Ad,h)
% K1 = J1(h) ;
