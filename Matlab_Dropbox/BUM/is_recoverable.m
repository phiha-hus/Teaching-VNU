function is_recoverable = is_recoverable(N,SRF, dist, s,eps)
    is_recoverable = 1;
    fc = floor(N/(2*SRF));
    n= 2*fc+1;
    
    Fn = dftmtx(N);
    Fn1 = Fn(1:(fc+1), :);
    Fn2 = Fn((N-fc+1):N, :);
    Fn = [Fn2;Fn1];
    
    support = find_support(N,SRF,dist,s);
    
    % test for a number signals
    for k = 1:20
        % construct a random sparse signal supported in 'support'
        x = zeros(N,1);
        x(support) = exp(1i* rand(1,length(support)));
        %x(support) = (-1).^randint(length(support),1) .* (10.^rand(length(support),1));
    
        % solve the basis pursuit to obtain an estimated signal
        cvx_begin
            variable z(N)
            Fn*z == Fn*x;
            minimize norm(z,1)
        cvx_end
        
        % if the error of the estimatate trepasses a toleranz, it means the
        % recovery is not successful
        if norm(z-x,1) > eps
            is_recoverable = 0;
        end
    end
end


