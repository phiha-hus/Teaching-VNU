function dydt = vidu23_fun(t,y,Z)
%DDEX1DE  Example of delay differential equations for solving with DDE23.
%
%   See also DDE23.

%   Jacek Kierzenka, Lawrence F. Shampine and Skip Thompson
%   Copyright 1984-2014 The MathWorks, Inc.

A     =  [-1 -3;2 -5]; 
Ad    =  [1.66 -0.697;0.93 -0.330];

dydt = A * y + Ad * Z;
