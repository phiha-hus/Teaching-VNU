function [value] = check_C2ctrb(M,D,K,B)
    [n,~]= size(M)
    value = 1;

    if(rank(M)~=n)
        value = 0;
        disp('this program cant test');
    else
        I = eye(n)
        o = zeros(n,n)
        N = inv(M)
        A =[[o,I]; [(-N*K) , (-N*D)]]
        nB =[0*B; N*B]
        if(rank([A,nB]) < n)
            value = 0
    end
end
