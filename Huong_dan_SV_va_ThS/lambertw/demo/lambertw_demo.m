% Demo script for lambertw

% Pascal Getreuer 2006

fprintf('\nIterated Exponentiation:\n');
fprintf(['\nIf the iterated exponentiation z^z^z^z^... converges, the limit\n',...
      'value can be expressed in term of the Lambert W function:\n',...
      '    z^z^z^z^... = -lambertw(-log(z))/log(z).\n',...
      'If z = 1.3, then the limit is computed to be\n']);

z = 1.3;
zlim = -lambertw(-log(z))/log(z)

fprintf(['To verify this numerically, z^z^z^z^...^z is computed out to\n',...
      '40 iterations.  The exponentiated value converges to zlim:\n']);
   
zz(1) = z;
zlim = -lambertw(-log(z))/log(z);

fprintf('  k    zz     zlim-zz\n');

for k = 1:40
   zz(k+1) = z^zz(k);
   
   if k < 10 | ~rem(k,10)
      fprintf(' %2d  %.4f   %.2g\n',k,zz(k+1),zlim-zz(k+1));
   end   
end

fprintf('\n');

figure;
plot(zz,'.-');
xlabel('Iteration');
ylabel('z\^z\^z\^...\^z');
xlim([1,25]);
title(sprintf('Iterated exponentiation of z = %g',z));

fprintf('\nW_0(z) Surface Plot:\n');
fprintf(['\nW_0(z) is the principal branch of the Lambert W function.  Its\n',...
      'magnitude in the complex plane is displayed in the surface plot.\n',...
      'The key features to observe is the point W_0(0) = 0 and the\n',...
      'roughly radial symmetry.\n\n']);

figure;
N = 73;
[x,y] = meshgrid(linspace(-6,6,N),linspace(-6,6,N));
w = abs(lambertw(x + 1i*y));

surf(x,y,w);
shading interp
colormap(hot(256));
hold on

n1 = 4;
n2 = 4;

for k = 1:n1:N
   plot3(x(k,:),y(k,:),w(k,:),'k-');
end

for k = 1:n2:N
   plot3(x(:,k),y(:,k),w(:,k),'k-');
end

axis([-6,6,-6,6,0,2.5]);
view(40,32);
xlabel('Re z');
ylabel('Im z');
title('|W_0(z)|');
set(gca,'XTick',[-6,0,6],'YTick',[-6,0,6],'ZTick',[0,2.5]);


