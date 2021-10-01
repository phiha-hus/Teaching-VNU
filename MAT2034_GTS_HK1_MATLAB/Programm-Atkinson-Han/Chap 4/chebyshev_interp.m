function [nodes, fcn_values, div_diff_fcn] = chebyshev_interp(n)

% This creates an interpolant of order n to the function
% fcn(x) on [-1,1],  which is given below as a function 
% subprogram. The nodes are the Chebyshev zeroes of the 
% degree n+1 Chebyshev polynomial on [-1,1]. The program 
% gives two plots: first the true function and its 
% interpolant, and second, the error in the interpolation.

% Create the nodes and associated divided differences.
h = pi/(2*(n+1));
nodes = cos(h*[1:2:2*n+1]);
fcn_values = fcn(nodes);
div_diff_fcn = divdif(nodes,fcn_values);

% Create the points at which the functions are to be 
% graphed.
x_eval = -1:.002:1;
true_fcn = fcn(x_eval);
y_eval = interp(nodes,div_diff_fcn,x_eval);

% Create the graph of the function and its interpolant.
m = min([min(true_fcn),min(y_eval)]);
if m > 0
    m = .9*m;
else
    m = 1.1*m;
end
M = max([max(true_fcn),max(y_eval)]);
if M < 0
    M = .9*M;
else
    M = 1.1*M;
end
axis([-1.1,1.1,m,M])
hold on
plot(x_eval,true_fcn,'r','LineWidth',1)
plot(x_eval,y_eval,':')
legend('True function','Interpolant')
hold on
plot(nodes,fcn_values,'.','MarkerSize',6)
hold off

pause
clf

% Create the window for the graph of the error.
error = true_fcn - y_eval;
M = max(error);
if M < 0
    M = .9*M;
else
    M = 1.1*M;
end
m = min(error);
if m > 0
    m = .9*m;
else
    m = 1.1*m;
end
axis([-1.1,1.1,m,M])
hold on

% Create the graph of the error in the interpolant.
plot(x_eval,error,'r','LineWidth',1)
hold off

% Print the maximum error.
disp(['maximum error = ',num2str(max(abs(error)))])

function fval = fcn(x)

fval = exp(x);
