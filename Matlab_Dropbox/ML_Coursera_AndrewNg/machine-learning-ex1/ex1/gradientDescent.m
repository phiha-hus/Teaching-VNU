function [theta, J_history] = gradientDescent(X, y, theta, alpha, num_iters)
%GRADIENTDESCENT Performs gradient descent to learn theta
%   theta = GRADIENTDESCENT(X, y, theta, alpha, num_iters) updates theta by 
%   taking num_iters gradient steps with learning rate alpha

% Initialize some useful values
m = length(y); % number of training examples
J_history = zeros(num_iters, 1);

    %theta = [-1;2]; % Initialize theta by Phi  
    %theta = [0; 0];
    
for iter = 1:num_iters

    % ====================== YOUR CODE HERE ======================
    % Instructions: Perform a single gradient step on the parameter vector
    %               theta. 
    %
    % Hint: While debugging, it can be useful to print out the values
    %       of the cost function (computeCost) and gradient here.
    %
    h_theta = X * theta;  
    [~,n] = size(X);
    %for j = 1:n
    %  theta(j) = theta(j) - alpha/m * dot(h_theta-y,X(:,j));
    %end
    
    theta = theta - alpha/m * (X'*(h_theta-y));

    % ================================  ============================

    % Save the cost J in every iteration    
    J_history(iter) = computeCost(X, y, theta);

end

end
