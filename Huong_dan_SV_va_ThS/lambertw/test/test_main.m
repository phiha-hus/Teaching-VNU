%The test_main test script for lambertw

% Pascal Getreuer 2006

fprintf('\n *** Inversion Error Test ***\n\n');
fprintf(['Ideally, lambertw(b,z)*exp(lambertw(b,z)) = z for any complex z\n',...
      'and any integer branch index b, but this is limited by machine\n',...
      'precision.  The inversion error |lambertw(b,z)*exp(lambertw(b,z)) - z|\n',...
      'is small but worth minding.\n\n',...      
      'Experimentation finds that the error is usually on the order of \n',...
      '|z|*1e-16 on the principal branch.  This test computes the inversion\n',...
      'error over the square [-10,10]x[-10,10] in the complex plane, large\n',...
      'enough to characterize the error away from the branch points at\n',...
      'z = 0 and -1/e.\n\n']);

N = 81;    % Use NxN points to sample the complex plane
R = 10;    % Sample in the square [-R,R]x[-R,R]
x = linspace(-R,R,N);
y = linspace(-R,R,N);
[xx,yy] = meshgrid(x,y);
z = xx + 1i*yy;

for b = -4:4
   w = lambertw(b,z);
   InvError = abs(w.*exp(w) - z);   
   fprintf('Largest error for b = %2d:  %.2e\n',b,max(InvError(:)));
   
   if b == 0
      figure;
      imagesc(x,y,InvError)
      axis image
      colorbar
      colormap(gray)
      title('Inversion error for W_0(z)');
      xlabel('Re z');
      ylabel('Im z');      
   elseif b == -1
      figure;
      imagesc(x,y,InvError);
      axis image      
      colorbar      
      colormap(gray)
      title('Inversion error for W_-_1(z)');
      xlabel('Re z');
      ylabel('Im z');
   end
end

fprintf('\n');
