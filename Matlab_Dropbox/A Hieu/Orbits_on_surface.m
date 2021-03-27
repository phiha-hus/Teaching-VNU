clear all
close all

% FIRST, lets plot the surface of optimal solution : 
% You need to run the script "SURFACE_.mat" before
load SURFACE.mat

% SURFACE e_3 : THE SURFACE TO WHICH THE ORBITS CONVERGE WITH OPTIMAL COSTS
fill3(points(1,:),points(2,:),points(3,:),'b')
grid on
alpha(0.3)
hold on


% % SURFACE e_3_rk4 : IF YOU WANT TO PLOT THIS SURFACE REMOV THE COMMENTS
% fill3(points_rk4(1,:),points_rk4(2,:),points_rk4(3,:),'r')
% grid on
% alpha(0.3)


%% Identifions les param : 
[K1,K2,K3,S1,S2,S3,c1,c2,c3,a,r,q,phi,p,epsilon] = getparam_3pays(1);

% !! PLEASE CHECK THE VALUE OF c1, c2, c3 in getparam_3pays.m

%% SIMU

t_fin = 100/epsilon; %temps rapide (= année/epsilon)



Nb_orbits_to_plot = 15
for i = 1:Nb_orbits_to_plot


% RANDOM INITIAL CONDITIONS : CI = [B1_init B2_init B3_init e1_init e2_init e3_init]    
    CI = [K1*rand(1) K2*rand(1) K3*rand(1) 100*S1*rand(1) 100*S2*rand(1) 100*S3*rand(1) ];    

% SOLVE THE DIFFERENTIAL EQUATIONS :     
tol = 1e-9*[1 1 1 1 1 1];
options = odeset('RelTol',1e-9,'AbsTol',tol);
[T,Y] = ode45(@rigid,[0 t_fin],CI,options);

B1 = Y(:,1);
B2 = Y(:,2);
B3 = Y(:,3);

e1 = Y(:,4);
e2 = Y(:,5);
e3 = Y(:,6);


% PLOT THE 3D ORBIT CORRESPONDING WITH THESE INITIAL CONDITIONS: 
a = lines(Nb_orbits_to_plot);
col = a(i,:);
figure(1)
plot3(e1, e2, e3,'color',col)
hold on
% Put a circle at final position: 
plot3(e1(end), e2(end),e3(end),'o','color',col,'linewidth',3)
% Put a cross at initial position: 
%plot3(e1(1), e2(1),e3(1),'+','color',col,'linewidth',3)

grid on


end

xlabel('Fishing Effort Zone 1','fontsize',18)
ylabel('Fishing Effort Zone 2','fontsize',18)
zlabel('Fishing Effort Zone 3','fontsize',18)
set(gca,'fontsize',15)