N = 4096; 
SRF = 8;
fc = floor(N/(2*SRF));
n = 2*fc+1;
eps = 1e-2;
s = 2; % size of T

T = zeros(s,0); 

% binary search for minimum distance
dist1 = SRF/2;
dist2 = 4*SRF;

while (dist2-dist1) > 2  
    dist = floor((dist1+dist2)/2);
    if is_recoverable(N,SRF,dist,s,eps) == 0
        dist1 = floor((dist1+dist2)/2)
    else
        dist2 = floor((dist1+dist2)/2)
    end
end



