% dist = minimum distance between two consecutives spikes
% s = number of spikes

function support = find_support(N,SRF,dist,s)
    % the given setting
    fc = floor(N/(2*SRF));
    n = 2*fc+1;
    
    Fn = dftmtx(N);
    Fn = Fn([(N-fc+1):N, 1:(fc+1)], :);    
    
    % initial value of support
    % initial value of the auxiliary matrix which contains all
    % selected columns of Fn
    support = zeros(1, s);
    A = zeros(n,0);
    
    % available indices at each turn
    a_index = 1:N; 
    
    % Greedy algorithm: construct adversarial support with elements
    % separated by 'dist' by sequentially adding elements to the support.
    % Each new element is chosen to maximize the condition number of the
    % submatrix 'A' formed by the columns corresponding to the selected 
    % elements. 
    
    for k = 1:s
        % find max/argmax cond([A, Fn(:, a_index(j))]) -> max_value/ ind
        max_value = cond([A,Fn(:, a_index(1))]); 
        ind = 1; 
        for j = 2:length(a_index)
            if cond([A, Fn(:, a_index(j))]) > max_value
                ind = j;
                max_value = cond([A, Fn(:, a_index(j))]);
            end
        end
        
        % add index of the right column of Fn in support
        % this is the location of the next spike
        support(k) = a_index(ind);
        
        % update matrix 'A'
        A = [A, Fn(:, a_index(ind))];
        
        % update a_index
        index_to_del = (a_index(ind)-dist+1):(a_index(ind)+dist-1);
        a1 = index_to_del(index_to_del <= 0) +N;
        a2 = index_to_del(index_to_del > N) -N;
        a3 = index_to_del(index_to_del>0 & index_to_del<= N);
        to_del = [a1,a2,a3];
        a_index = setdiff(a_index, to_del);
    end
end        