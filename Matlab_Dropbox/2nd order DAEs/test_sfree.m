% The Test Set for solver sfree_triple and its Supporting Functions
% 1st Example: Example 3.1 , Lena Wunderlich's Dissertation
% 2nd Example: Modify Example 3.1, where more redundant equations are added
% 3rd Example: A simple 2 by 2 systems, to illustrate the difference of
% counting methods between Shi's approach and my approach.
% 4th Example: System of size 3 by 3
% 5th Example: System of size 3 by 3, but not regular (i.e. non-unique solution)
% 6th Example with random coefficients systems

function test_sfree(index)

switch index
  case 1
        %% Example 3.1, Lena Wunderlich's Dissertation
        A2 = [1 0;0 0]
        A1 = [1 0;0 0]
        A0 = [0 1;1 0]    
  case 2
        %% 2nd Example - Over-determined system
        A2 = [1 0;0 0;0 0;0 0]
        A1 = [0 0;1 0;0 0;0 0]
        A0 = [0 1;0 0;1 0;0 0]       
  case 3
        %% 3rd Example
        A2 = [1 0; 0 0];
        A1 = [0 0; 0 0];
        A0 = [0 0; 0 1];    
  case 4
        %% 4th Example: system of size 3 by 3
        A2 = [0 1 0; 0 0 1; 0 0 0]
        A1 = [0 0 0; 0 1 0; 0 0 1]
        A0 = [1 0 0; 0 0 0; 0 0 0]  
        h = 0.01
        A0 = h.^2 * A0 - h * A1 + A2
        A1 = h * A1 - 2 * A2
  case 5
        %% 5th Example: another system of size 3 by 3, but non-unique solution
        A2 = [0 1 0; 0 0 1; 0 0 0]
        A1 = [0 0 0; 0 1 0; 0 0 0]
        A0 = [0 0 0; 0 0 0; 0 0 1]  
  case 6
        %% 6th Example with random coefficients systems
        m = 10; n = 10;
        rm = 4; rn = 5; rp = 7;
        P = rand(m,m); Q = rand(n,n);
        A2 = P * [rand(rm,n); zeros(m-rm,n)] * Q;
        A1 = P * [rand(rn,n); zeros(m-rn,n)] * Q;
        A0 = P * [rand(rp,n); zeros(m-rp,n)] * Q;
    case 7     
        disp('Anh Numberone dep trai de thuong')  
end
  
%[is_sfree,M,N,P] = check_sfree_v1(A2,A1,A0)
%[is_sfree,M,N,P] = check_sfree_v2(A2,A1,A0)
[ka,M,N,P] = sfree_triple(A2,A1,A0)

