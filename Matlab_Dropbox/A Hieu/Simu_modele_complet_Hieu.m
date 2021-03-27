clear all
close all

%% Identifions les param : 
[K1,K2,K3,S1,S2,S3,c1,c2,c3,a,r,q,phi,p,epsilon] = getparam_3pays(1);


% densite de carcapa : 
k1 = K1/S1;
k2 = K2/S2;
k3 = K3/S3;

% Suivant ces param, quel est le cout max de l'effort de p?che dans chaque
% pays qui permet l'existance d'une p?cherie p?renne?
c1_max = p*q*k1;
c2_max = p*q*k2;
c3_max = p*q*k3;

% Suivant ces param, quel est le cout  de l'effort de p?che dans chaque
% pays qui permet d'optimiser les captures totales?
c1_opt = p*q*k1/2;
c2_opt = p*q*k2/2;
c3_opt = p*q*k3/2;


% CALCULONS LES B_barre : 
K = K1 + K2 + K3;
B1_barre = (K*c1*S1)/(p*q*K1);
B2_barre = (K*c2*S2)/(p*q*K2);
B3_barre = (K*c3*S3)/(p*q*K3);

% L'effort prevu a l'equilibre dans la zone 3 (car B3_barre est le plus
% petit) : 

E_eq3 = (r*S3*(K-B3_barre))/(q*K3)

% et donc Landings prevu a l'equilibre = 
E_eq3*q*B3_barre

%c1_opt/k1

%Pecherie qui survie : 
lite_ratio = [c1/k1 c2/k2 c3/k3];
find(lite_ratio == min(lite_ratio));
%% SIMU

t_fin = 100/epsilon; %temps rapide (= ann?e/epsilon)




%%
% CONDITIONS INITIALES : 
% Peche partout : 
CI = [10 1 1 100 0.11 0.01 ];

CI_rec=[];
Eq_rec = [];

% TEST SENSIB CONDITIONS INITIALES : 
CI_rec=[];
Eq_rec = [];

lim_E1 = 100;
lim_E2 = 100;
lim_E3 = 100;

nbline = 10
for i = 1:1%nbline%00
  %  close all
% Biom init entre 0 et K
% Effort init entre 0 et 365*S
%CI = [K1*rand(1) K2*rand(1) K3*rand(1) 100*S1*rand(1) 100*S2*rand(1) 100*S3*rand(1) ];    
%CI = [10*rand(1) 10*rand(1) 10*rand(1) 100*rand(1) 100*rand(1) 100*rand(1) ];    

CI = [ 1.8      1.2     0.9 90     80     74]
CI2 = [ 1.9      1.1     0.6 80     85     70]
CI3 = [ 1.7      0.8     0.9 70     95     60]

% Peche Maroc REK : 
%CI = [4 3 2 1 0 0];

tol = 1e-9*[1 1 1 1 1 1];

options = odeset('RelTol',1e-9,'AbsTol',tol);
[T,Y] = ode45(@rigid_new,[0 t_fin],CI,options);
[T2,Y2] = ode45(@rigid_new,[0 t_fin],CI2,options);
[T3,Y3] = ode45(@rigid_new,[0 t_fin],CI3,options);

B1 = Y(:,1);
B2 = Y(:,2);
B3 = Y(:,3);

e1 = Y(:,4);
e2 = Y(:,5);
e3 = Y(:,6);

t = T * epsilon;
figure(10); clf;
subplot(2,1,1)
plot(t,B1,'b',t,B2,'r',t,B3,'--k','linewidth',1.5)
ylabel('Fishing Biomass','fontsize',18)
xlabel('Years','fontsize',18)
xlim([0,100])
ylim([0,3])
grid on

subplot(2,1,2)
plot(t,e1,'b',t,e2,'r',t,e3,'--k','linewidth',1.5)
ylabel('Fishing Effort','fontsize',18)
xlabel('Years','fontsize',18)
xlim([0,100])
%ylim([60,150])
grid on

%print -depsc Fig3_3.eps

serie_temp = 0;
if serie_temp == 1
% Landings: 
L1 = q*B1.*e1./S1;
L2 = q*B2.*e2./S2;
L3 = q*B3.*e3./S3;

figure
temps = T*epsilon; %  TEMPS LENT

subplot(2,1,1)
plot(temps, B1,'b')
hold on
plot(temps, B2,'r')
plot(temps, B3,'k')

% plot(temps, L1,'--b')
% hold on
% plot(temps, L2,'--r')
% plot(temps, L3,'--k')
%axis([0 t_fin 0 max(n)+max(n)/10])
%legend('Morocco', 'Mauritania','Senegal')
legend('Zone 1', 'Zone 2','Zone 3')%, 'Landings zone 1', 'Landings zone 2', 'Landings zone 3')
ylabel('Fish Biomass')
grid on
title(['c_1 = ' num2str(c1),'    c_2 = ' num2str(c2),'    c_3 = ' num2str(c3)])

subplot(2,1,2)
plot(temps, e1,'b')
hold on
plot(temps, e2,'r')
plot(temps, e3,'k')
ylabel('Fishing Effort')
grid on
xlabel('Year')
title(['Effort Costs : [c_1 c_2 c_3 ] = ' num2str([c1 c2 c3]),' ;  Conditions Initiales : [B1 B2 B3] = ' num2str(CI(1:3)) '  ;  [E1 E2 E3] = ' num2str(CI(4:6))])

pos = [3         507        1438         299];
set(gcf,'position', pos)

end

a = lines(nbline);
col = a(i,:);
figure(1)
plot3(e1, e2, e3,'color',col)
hold on
plot3(e1(end), e2(end),e3(end),'o','color',col,'linewidth',3)
grid on
xlabel('Zone 1')
ylabel('Zone 2')
zlabel('Zone 3')

%print -depsc Fig1_3.eps


FIG2 = 0
if FIG2 == 1
figure(2)
subplot(3,1,1)
plot(B1, e1,'color',col)
hold on
plot(B1(end), e1(end),'o','color',col,'linewidth',3)
ylabel('Fishing Effort')
xlabel('Fish Biomass')
title('ZONE 1')
lim_E1 = max(lim_E1,max(e1));
axis([0 4 0 lim_E1])

subplot(3,1,2)
plot(B2, e2,'color',col)
hold on
plot(B2(end), e2(end),'o','color',col,'linewidth',3)
ylabel('Fishing Effort')
xlabel('Fish Biomass')
title('ZONE 2')
lim_E2 = max(lim_E2,max(e2));
axis([0 3 0 lim_E2])

subplot(3,1,3)
plot(B3, e3,'color',col)
hold on
plot(B3(end), e3(end),'o','color',col,'linewidth',3)
ylabel('Fishing Effort')
xlabel('Fish Biomass')
title('ZONE 3')
%grid on
lim_E3 = max(lim_E3,max(e3));
axis([0 2.2 0 lim_E3])

end

%Figure avec les 3 equilibres sur le m?me plot
FIG3=1
if FIG3 ==1
B = B1+B2+B3;
e = e1+e2+e3;
Bnew = Y2(:,1)+Y2(:,2)+Y2(:,3);
enew = Y2(:,4)+Y2(:,5)+Y2(:,6);
Bnew2 = Y3(:,1)+Y3(:,2)+Y3(:,3);
enew2 = Y3(:,4)+Y3(:,5)+Y3(:,6);

figure
% % EN GRIS : 
% plot(B1, e1,'k')
% hold on
% plot(B2, e2,'color',[0.5 0.5 0.5])
% plot(B3, e3,'k--')
% 
% %plot(B, e1+e2+e3,'--b')
% plot(B1(end), e1(end),'ok','linewidth',3)
% plot(B2(end), e2(end),'o','color',[0.5 0.5 0.5],'linewidth',3)
% plot(B3(end), e3(end),'ok','linewidth',3)
% 
% plot(B1(1), e1(1),'+k','linewidth',3)
% plot(B2(1), e2(1),'+','color',[0.5 0.5 0.5],'linewidth',3)
% plot(B3(1), e3(1),'+k','linewidth',3)

% EN COULEUR : 
plot(B1, e1,'b'); grid on
hold on
plot(B2, e2,'r'); grid on
plot(B3, e3,'--k'); grid on

%plot(B, e1+e2+e3,'--b')
plot(B1(end), e1(end),'ob','linewidth',2)
plot(B2(end), e2(end),'or','linewidth',2)
plot(B3(end), e3(end),'ok','linewidth',2)

plot(B1(1), e1(1),'+b','linewidth',2) 
plot(B2(1), e2(1),'+r','linewidth',2)
plot(B3(1), e3(1),'+k','linewidth',2)

ylabel('Fishing Effort','fontsize',18)
xlabel('Fish Biomass','fontsize',18)
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.7;
%ax.Fontsize = 10;
%ax.Linewidth = 0.8;
grid on

print -depsc Fig2_3.eps

%legend('Zone 1', 'Zone 2','Zone 3')%, 'TOTAL')
%title(['CFUE: [c_1 c_2 c_3 ] = ' num2str([c1 c2 c3]),' ;  Conditions Initiales : [B1 B2 B3] = ' num2str(CI(1:3)) '  ;  [E1 E2 E3] = ' num2str(CI(4:6))])
%title(['CFUE: c_1  = ' num2str(c1) ' ; c_2 = '  num2str(c2) ' ; c_3 = '  num2str(c3) ' ;  B1_{INIT}  = ' num2str(CI(1)) ' ; B2_{INIT} = ' num2str(CI(2)) ' ; B3_{INIT} = ' num2str(CI(3)) '  ;  E1_{INIT} = ' num2str(CI(4))  ' ; E2_{INIT} = ' num2str(CI(5)) ' ; E3_{INIT}  = ' num2str(CI(6))])
%title(['CFUE: c_1  = ' num2str(c1) ' ; c_2 = '  num2str(c2) ' ; c_3 = '  num2str(c3) ])
lim_E = max([max(e1) max(e2) max(e3)]);
lim_B = max([max(B1) max(B2) max(B3)]);
set(gca,'fontsize',16)

%axis([0 lim_B 0 lim_E])
display(['t = 100 : [B1 B2 B3] = ' num2str([B1(end) B2(end) B3(end)])  ' [E1 E2 E3] = ' num2str([e1(end) e2(end) e3(end)])]);

% pos = [ 134   235   697   571];
% set(gcf,'position', pos)
end

% Catch
Catch1 = B1(end)*e1(end)*q;
Catch2 = B2(end)*e2(end)*q;
Catch3 = B3(end)*e3(end)*q;

display(['Total Catch = ' num2str(Catch1)  ' + ' num2str(Catch2)  '+ ' num2str(Catch3)  ' = ' num2str(Catch1+Catch2+Catch3) ]);

CI_rec = [CI_rec; CI];
Eq_rec = [ Eq_rec ; [B1(end) B2(end) B3(end) e1(end) e2(end) e3(end)] ];


figure(4); clf;
plot(B,e,'-k')
hold on

plot(Bnew,enew,'-k'); hold on 

plot(Bnew2,enew2,'-k'); hold off
%legend('1st Initial Condition','2nd Initial Condition','3rd Initial Condition')
ylabel('Total Fishing Effort','fontsize',18)
xlabel('Total Fish Biomass','fontsize',18)
grid on

%print -depsc Fig4_q.eps

end


%%
pause

figure
load SURFACE_rksur4.mat
fill3(points(1,:),points(2,:),points(3,:),'r')
grid on
alpha(0.3)
hold on
load SURFACE.mat

fill3(points(1,:),points(2,:),points(3,:),'b')
grid on
alpha(0.3)


figure
liste_eq_B1 = Eq_rec(:,1);
liste_eq_B2 = Eq_rec(:,2);
liste_eq_B3 = Eq_rec(:,3);

liste_eq_e1 = Eq_rec(:,4);
liste_eq_e2 = Eq_rec(:,5);
liste_eq_e3 = Eq_rec(:,6);

plot(liste_eq_B1,liste_eq_e1,'b.', 'markersize',20)
hold
plot(liste_eq_B2,liste_eq_e2,'r.', 'markersize',20)
plot(liste_eq_B3,liste_eq_e3,'k.', 'markersize',20)

ylabel('Fishing Effort','fontsize',21)
xlabel('Fish Biomass','fontsize',21)
grid on
legend('Zone 1', 'Zone 2','Zone 3')%, 'TOTAL')
title(['100 Equilibriums according to the initial conditions, with optimal costs'],'fontsize',17)
set(gca,'fontsize',19)



%% les catch

figure
for i = 1:3
    Catch(:, i) = Eq_rec(:,i).*Eq_rec(:,i+3).*q;
end

bar(Catch,'stack')

xlabel('Simulation')
ylabel('Landings')
legend('Zone 1', 'Zone 2','Zone 3')%, 'TOTAL')
grid on

%% EFFET DE L'EFFORT INITIAL DE PECHE
figure
plot(CI_rec(:,4), Eq_rec(:,4),'.','markersize',20)
hold
plot(CI_rec(:,5), Eq_rec(:,5),'r.','markersize',20)
plot(CI_rec(:,6), Eq_rec(:,6),'k.','markersize',20)
xlabel('Initial Fishing Effort','fontsize',21)
ylabel('Fishing Effort at Equilibrium','fontsize',21)
legend('Zone 1', 'Zone 2','Zone 3')%, 'TOTAL')

title('Impact of Initial Fishing Effort on the Equilibrium with Optimal costs','fontsize',17)
grid on
set(gca,'fontsize',19)


%% SERIE TEMPORELLE
figure
% % GRAY : 
% subplot(2,1,1)
% plot(T*epsilon,Y(:,1),'k')
% hold on
% plot(T*epsilon,Y(:,2),'color',[0.5 0.5 0.5])
% plot(T*epsilon,Y(:,3),'k--')
% 
% %legend('Zone 1', 'Zone 2','Zone 3')%, 'TOTAL')
% ylabel('Fish Biomass','fontsize',18)
% %xlabel('Years')
% grid on
% set(gca,'fontsize',16)
% 
% subplot(2,1,2)
% plot(T*epsilon,Y(:,4),'k')
% hold on
% plot(T*epsilon,Y(:,5),'color',[0.5 0.5 0.5])
% plot(T*epsilon,Y(:,6),'k--')
% %legend('Zone 1', 'Zone 2','Zone 3')%, 'TOTAL')
% ylabel('Fishing Effort','fontsize',18)
% xlabel('Years','fontsize',18)
% grid on
% set(gca,'fontsize',16)

% COULEUR : 
subplot(2,1,1)
plot(T*epsilon,Y(:,1),'b')
hold on
plot(T*epsilon,Y(:,2),'r')
plot(T*epsilon,Y(:,3),'g')

%legend('Zone 1', 'Zone 2','Zone 3')%, 'TOTAL')
ylabel('Fish Biomass','fontsize',18)
%xlabel('Years')
grid on
set(gca,'fontsize',16)

subplot(2,1,2)
plot(T*epsilon,Y(:,4),'b')
hold on
plot(T*epsilon,Y(:,5),'r')
plot(T*epsilon,Y(:,6),'g')
%legend('Zone 1', 'Zone 2','Zone 3')%, 'TOTAL')
ylabel('Fishing Effort','fontsize',18)
xlabel('Years','fontsize',18)
grid on
set(gca,'fontsize',16)
