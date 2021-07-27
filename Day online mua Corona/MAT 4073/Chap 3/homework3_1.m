clear; clc;

% Exercíse 3.5
disp('Exercise 3.5')
% spline(x,y,z)
% - x is some node
% - y = f(x)
% - z: output = f(z)
years = [0,5,10,15];

westEuro = [72.8, 74.2, 75.2, 76.4];
eastEuro = [70.2, 70.2, 70.3, 71.2];
disp('life expectation in 1977, 1983, 1988 of weast:')
spline(years, westEuro, [2,8,13])
disp('life expectation in 1977, 1983, 1988 of east:')
spline(years, eastEuro, [2,8,13])

% Exercise 3.6
disp('Exercise 3.6')
format shorte
t = [4, 8, 12, 16, 20];
p = [1000.7794, 1000.6427, 1000.2805, 999.7165, 998.9700];

tTest = [6, 10, 14, 18];
pTest = [1000.74088, 1000.4882, 1000.0224, 999.3650];

disp('Value in [6,10,14,18]')
spline(t,p,tTest)

disp('Error of p in [6, 10, 14, 18]:')
sum(abs(spline(t,p,tTest) - pTest))

% Exercise 3.7
disp('Exercise 3.7')
year = [1965, 1970, 1980, 1985, 1990, 1991];
product = [17769, 24001, 25961, 34336, 29036, 33417];
time = [1962, 1977, 1992];
realValues = [12380, 27403, 32059];

disp('Values in [1962, 1977, 1992]')
values = spline(year, product, time)
disp('Error values with real values:')
sum(abs(values - realValues))


