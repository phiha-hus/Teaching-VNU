% Test real part of the spectrum using Lambert functions

clear all; close all; clc

%%
% Test the eigenvalues correspond to the stability of DDAEs in the new
% sense (weakly exponentially stable)
N = 10; Spec = [];

gamma = 0.5;
% gamma = rand(1,1)
y = -1/gamma;

for i=-N:1:N
    % disp('The real part of the eigenvalues on the k-branch of the Lambert function w.r.t x=-1 is ')
    lambda = -real(lambertw(i,y));
    Spec = [Spec lambda];    
end

subplot(2,1,1)
plot(-N:1:N,Spec,'*r'); 
%title('Real Parts of Eigenvalues from branches -10 to 10 of Lambert W function')
xlabel('Branches of Eigenvalues'); 
ylabel('Real Parts');
grid on 

print -depsc Exam2_Ha18Acta_LambertW.eps
  ! epstopdf Exam2_Ha18Acta_LambertW.eps
  ! rm Exam2_Ha18Acta_LambertW.eps
 

%%
for gamma = -0.9:0.1:0.9
  for k=-N:N
    K = [K k];
    % disp('The real part of the eigenvalues on the k-branch of the Lambert function w.r.t x=-1 is ')
    y = -1/gamma;
    lambda = real(lambertw(k,y));
    if lambda>=0
        Spec = [Spec lambda];
    else
        gamma
        k
        lambda        
        return
    end
  end
end



