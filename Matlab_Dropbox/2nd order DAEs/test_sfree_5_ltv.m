% The Test Set for solver sfree_triple_LTV and its Supporting Functions

%% Example 3.1, Lena Wunderlich's Dissertation

clear all; close all; clc

A2 = @(n)[1 0;0 0];

A1 = @(n)[1 0;0 0];

A0 = @(n)[0 1;1 0];

[ka,M,N,P]=sfree_triple_LTV(1,A2,A1,A0);

%% Example 3.14, Lena Wunderlich's dissertation

clear all; close all; clc

A2 = @(n)[1 0 0 0; 0 n 0 0; 0 0 0 0; 0 0 0 0];

A1 = @(n)[0 0 0 n; n 0 0 0; 0 0 1 0; 0 0 0 0];

A0 = @(n)[0 0 0 0; 0 0 0 0; 1 0 0 0; 0 0 1 0];

[ka,M,N,P]=sfree_triple_LTV(1,A2,A1,A0);


%% Example 3.23, Lena Wunderlich's dissertation

clear all; close all; clc

A2 = @(n)[n 0 0; 0 1 1; 0 n n];

A1 = @(n)[1 0 0; 0 0 0; 0 0 0];

A0 = @(n)[1 0 0; 0 1 0; 0 1+n 1];

[ka,M,N,P]=sfree_triple_LTV(1,A2,A1,A0);


Example 3.36