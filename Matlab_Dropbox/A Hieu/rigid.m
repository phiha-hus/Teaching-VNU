function dy = rigid(t,y)
dy = zeros(6,1);    % a column vector

[K1,K2,K3,S1,S2,S3,c1,c2,c3,a,r,q,phi,p,epsilon] = getparam_3pays(t);

% Pour facilite l'ecriture, on associe les bon noms de variables aux y(1),
% y(2), ... :
B1 = y(1);
B2 = y(2);
B3 = y(3);

e1 = y(4);
e2 = y(5);
e3 = y(6);


% MOdele COMPLET:
dB1 = (a*B2)/K2 - (a*B1)/K1               + epsilon*( r*B1*(1-B1/K1)   - (q*B1*e1)/S1 );
dB2 = (a*B1)/K1 + (a*B3)/K3 - (2*a*B2)/K2 + epsilon*( r*B2*(1 - B2/K2) - (q*B2*e2)/S2 );
dB3 = (a*B2)/K2 - (a*B3)/K3               + epsilon*( r*B3*(1-B3/K3)   - (q*B3*e3)/S3 );

de1 = epsilon*(phi/c1)* (p*q*e1*(B1/S1) - c1*e1);
de2 = epsilon*(phi/c2)* (p*q*e2*(B2/S2) - c2*e2);
de3 = epsilon*(phi/c3)* (p*q*e3*(B3/S3) - c3*e3);


dy(1) = dB1; 
dy(2) = dB2; 
dy(3) = dB3;

dy(4) = de1; 
dy(5) = de2; 
dy(6) = de3;


