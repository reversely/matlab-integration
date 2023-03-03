%% MECH 222 Computer Lab 2

%% Function 1
%  f(x) = 1/(1 + x^2) over [0,1]
clear; close; clc;
a = 0; b = 1; N = 10; dx = (b - a)/N;
x = linspace(a,b,N + 1);
y = 1./(1 + x.^2);
plot(x,y,'.-')

% Right Riemann sum
Iright = sum(y(2:end))*dx

% Left Riemann sum
Ileft = sum(y(1:end-1))*dx

% Midpoint Riemann sum
xmid = linspace(a + dx/2,b - dx/2,N);
ymid = 1./(1 + xmid.^2);
Imid = sum(ymid)*dx

% Exact value
I = pi/4

%% Function 2
%  f(x) = sin(cos(x)) over [0,pi]

%% Find upper bound of |f'(x)| analytically
%  f'(x) = -cos(cos(x))*sin(x) => |f'(x)| <= 1

%% Find upper bound of|f'(x)| graphically
clear; close; clc;
X = linspace(0,pi,100);
Y = -cos(cos(X)).*sin(X);
plot(X,Y)

%% Approximate and estimate error
clear; close; clc;
K1 = 1; a = 0; b = pi/2; N = 100; dx = (b - a)/N;
x = linspace(a,b,N + 1);
y = sin(cos(x));
plot(x,y,'.-')

Iright = sum(y(2:end))*dx
error = (b - a)^2/2/N*K1
