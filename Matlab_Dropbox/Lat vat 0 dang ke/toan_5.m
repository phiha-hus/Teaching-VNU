% Test sum of the digits of a natural number
% When will it leave the set {1,4,7,9}

nmax = 100000;
Array1 =[1 2 3];
Array2 = [1 4 9];

for i=4:nmax
    s = i^2;    
    while s>=10
        s = sumOfDigit(s,10)
    end
    
    Array1 = [Array1 i];
    Array2 = [Array1 s];
    
    if (s~=1) && (s~=4) && (s~=7) && (s~=9) 
        disp('The result is')
        Array1
        Array2
        break
    end
end