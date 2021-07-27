for n = 4:13
    n 
    tic
    A = hilb(n)
    cond(A)
    U = chol(A)
    [Q,R] = qr(A)
    [U,S,V] = svd(A)
    toc
end
