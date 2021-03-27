clear all
close all

%% PARAM % Identifions les param : 
[K1,K2,K3,S1,S2,S3,c1,c2,c3,a,r,q,phi,p,epsilon] = getparam_3pays(1);

% CALCULONS LES B_barre : 
K = K1 + K2 + K3;
B1_barre = (K*c1*S1)/(p*q*K1);
B2_barre = (K*c2*S2)/(p*q*K2);
B3_barre = (K*c3*S3)/(p*q*K3);


% LES ALPHA
alpha1 = K1/K;
alpha2 = K2/K;
alpha3 = K3/K;

k1 = K1/S1;
k2 = K2/S2;
k3 = K3/S3;

k = k1+k2+k3;
for e_1 = 0:500;
    for e_2 = 0:500;
e_1
       e_3(e_1+1,e_2+1) = (S3/(q*alpha3)) * ( r*(1-B1_barre/K) - ((q*alpha1)/S1 )*e_1 - ((q*alpha2)/S2 )*e_2 );
       e_3_rk4(e_1+1,e_2+1) = (S3/(q*alpha3)) * ( (r*K)/4 - ((q*alpha1)/S1 )*e_1 - ((q*alpha2)/S2 )*e_2 );
    end    
end

e_3(find(e_3<0))=nan;
e_3_rk4(find(e_3_rk4<0))=nan;

% figure
% surf(0:500,0:500,e_3)
% 
% xlabel('e1')
% ylabel('e2')
% zlabel('e3')


figure

% SURFACE e_3 : 
pointA = [0 0 max(max(e_3))]
[I J] = find(e_3>0);
pointB = [max(I) 0 0]
pointC = [0 max(J) 0]
points=[pointA' pointB' pointC']; 
fill3(points(1,:),points(2,:),points(3,:),'b')
grid on
alpha(0.3)
hold on

% SURFACE e_3_rk4 : 
pointA = [0 0 max(max(e_3_rk4))]
[I J] = find(e_3_rk4>0);
pointB = [max(I) 0 0]
pointC = [0 max(J) 0]
points_rk4=[pointA' pointB' pointC']; 
fill3(points_rk4(1,:),points_rk4(2,:),points_rk4(3,:),'r')
grid on
alpha(0.3)


xlabel('Fishing Effort Zone 1','fontsize',18)
ylabel('Fishing Effort Zone 2','fontsize',18)
zlabel('Fishing Effort Zone 3','fontsize',18)

set(gca,'fontsize',15)

save SURFACE
