figure(1); clf;

subplot(2,2,1)
%semilogy(T,sol_exact(1,:),'r -.')
plot(T,sol_exact(1,:),'r -.')
title('Exact solution 1st term by dde23')
grid on 

subplot(2,2,2)
%semilogy(T,sol_exact(2,:),'b -.')
plot(T,sol_exact(2,:),'b -.')
title('Exact solution 2nd term by dde23')
grid on

subplot(2,2,3)
%semilogy(T,sol_approx(1,:),'r -.')
plot(T,sol_approx(1,:),'r -.')
title('Approx solution 1st term')
grid on 

subplot(2,2,4)
%semilogy(T,sol_approx(2,:),'b -.')
plot(T,sol_approx(2,:),'b -.')
title('Approx solution 2nd term')
grid on


figure(2); clf;
subplot(2,2,1)
semilogy(T,err(1,:),'r -.')
title('Absolute error 1st term')
grid on 

subplot(2,2,2)
semilogy(T,err(2,:),'b -.')
title('Absolute error 2nd term')
grid on

subplot(2,2,3)
semilogy(T,err(1,:)./sol_exact(1,:),'r -.')
title('Relative error 1st term')
grid on 

subplot(2,2,4)
semilogy(T,err(2,:)./sol_exact(2,:),'b -.')
title('Relative error 2nd term')
grid on


figure(3); clf;
subplot(2,2,1)
plot(T,sol_approx(1,:),'r -','Linewidth',2); hold on
plot(T,sol_approx(2,:),'b -.','Linewidth',2); hold off
legend('x_1','x_2')
grid on 

subplot(2,2,2)
plot(T,err(1,:),'r -','Linewidth',2); hold on 
plot(T,err(2,:),'b -.','Linewidth',2); hold off
legend('Abs. err. x_1','Abs. err. x_2')
grid on 

print(gcf, '-dpdf', 'hinh3.pdf');