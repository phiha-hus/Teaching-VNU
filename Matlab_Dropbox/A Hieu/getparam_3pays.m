function [K1,K2,K3,S1,S2,S3,c1,c2,c3,a,r,q,phi,p,epsilon] = getparam_3pays(t)

% MODELE COMPLET PECHE 3 PAYS HIEU
% % PARAM IDENTIFIE : 
K1 = 4.3; %1e6 tonnes
K2 = 3.3; %1e6 tonnes
K3 = 2.4; %1e6 tonnes

% Surfaces des zones (1e4 km2)
S1 = 4.1472; 
S2 = 2.2720;
S3 = 1.5040;

a = 1; % megatonnes par 10 jours qui transitent d'une zone ? l'autre (au total K = 10 megatonnes)

r = 0.88;

q = 0.01;
% essais dec2017 :
%q = 0.02;%1e-3%1e-4%0.01; % "fraction of the average biomass in given area catch by a fishing vessel trawling on the entire area" No DIMMENSION

phi = 1; % taux de reinvestissement par ans. autour de 1.

% Prix du poisson (millions $ /megatone)
%p = 500 ;
%epsilon = 1/365; %; 10/365;

% essais dec2017 :
p = 500; % millions USD per 10000 km2 fished/prospected
epsilon = 10/365;




%% cout d'effort de p?che (par 1e4 km2 par ans) :
f = 1; % millons de dollars par 1000m3 de fuel)

% densite de carcapa : 
k1 = K1/S1;
k2 = K2/S2;
k3 = K3/S3;


% Compute max cost for the surviving of the 3 fisheries : 
% Suivant ces param, quel est le cout max de l'effort de p?che dans chaque
% pays qui permet l'existance d'une p?cherie p?renne?
c1_max = p*q*k1;
c2_max = p*q*k2;
c3_max = p*q*k3;

%Compute optimal cost for that maximise the total catch :
% Suivant ces param, quel est le cout  de l'effort de p?che dans chaque
% pays qui permet d'optimiser les captures totales?
c1_opt = p*q*k1/2;
c2_opt = p*q*k2/2;
c3_opt = p*q*k3/2;

% UNIFORM COST : only one winner
% c1 = 3*f;
% c2 = 3*f;
% c3 = 3*f;


% the order of magnitude of c can be increased if you take a smaller q... 

% essais dec2017 :
c1 = 3;
c2 = 3;
c3 = 3;

% OPTIMAL COST : THE 3 FISHERIES CO-EXIST (BUT NOT OPTIMAL HAVESTING???)
% c1 = c1_opt;
% c2 = c2_opt;
% c3 = c3_opt;


