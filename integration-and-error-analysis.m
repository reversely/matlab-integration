%  f(x) over [a,b]
clear; close; clc;

% Error needs to be less than (b-a)^3/(12N^2)*K2 < 10^-3.

% We want our k2 to be of the largest magnitude possible
% within the bounds of the equation |f"(x)| <= K2.
% Our function, f=@(x)sin(x.^2), has a double derivative of
% 2cos(x^2)-4x^2sin(x^2).
% Now when plotting the function from given range, 0 to sqrt(pi/2),

a = 0; b = sqrt(pi/2);

% the largest magnitude in is 2pi.

k2 = 2*pi;

% We find n by plugging in other known values into the error formula
%(b-a)^3/(12n)^2*k2 and solving for n.
% Plugging in values: n = (sqrt(pi/2)-(0)/(12*N^2)*2pi < 10^-3
% "ceil" rounds the value up to the highest integer.

n = ceil(sqrt((2*pi/10^-3)*(pi/2)^(3/2)*(1/12)))

% Now, we can use this value of n to return a trapezoid approximation function.
xk = linspace(a,b,n+1);
dx = (b-a)/n;

f=@(x)sin(x.^2);
fxk = sum(f(xk(2:n))); fxksub = sum(f(xk(1:n-1)));

trap_approx = dx*1/2*(fxk+fxksub)
error = (b-a)^3/(12*n^2)*k2
