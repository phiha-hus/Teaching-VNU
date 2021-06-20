clear all ; close all ; clc

%% Vi du 1.1

z = rand

k = -1:2

w = lambertw(k,z)


%%
A = [1 1 1;0 1 2;0 0 3];

[V,J] = jordan(A) ;

inv(V) * A * V - J ;

% WJ = [lambertw(0,3) 0 0 ; 0 lambertw(2,1) lambertw(2,1)/(1 + lambertw(2,1)) ; 0 0 lambertw(2,1)]

k = 2 ; 

WJ = [lambertw(k,3) 0 0 ; 0 lambertw(k,1) lambertw(k,1)/(1 + lambertw(k,1)) ; 0 0 lambertw(k,1)] ;

lambertw_matrix(k,J) ;

% WA = V * WJ * inv(V)





