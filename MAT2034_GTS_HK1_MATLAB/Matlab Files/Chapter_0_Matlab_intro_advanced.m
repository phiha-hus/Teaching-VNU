%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% More advanced topic
 
memory
error('Message for error: I donot know why ')

%% Machine precision

% Specify out precision in Octave
%output_precision             % Octave
%output_precision(7)          % Octave

% Specify out precision in Matlab
format short
format long
format shortE
format longE

%% Realmin, realmax, under and overfolows, floating points
realmin
realmax 
a = 1.0e+308;
b = 1.1e+308;
c = -1.01e+308;
(a+b)+c
a+ (b+c)

%% Lost of significant digits %See Quarteroni book

x = 1e-8
sin(x) - x
x.^3

%% Number of argument in/out
nargin   % Remark: nargin function can be used both internal and externally of a function. Check nargout for matrix/vector then see
nargout
