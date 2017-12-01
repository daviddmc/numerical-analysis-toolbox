function [x, y] = PredictorCorrector( f, a, b, y0, n )
%PREDICTORCORRECTOR Summary of this function goes here
%   Detailed explanation goes here

h = (b - a) / n;
x = zeros(n+1,1);
y = zeros(n+1,1);
x(1) = a;
y(1) = y0;

for i = 2:4
    k1 = h * f(x(i-1), y(i-1));
    k2 = h * f(x(i-1) + h/2, y(i-1) + k1/2);
    k3 = h * f(x(i-1) + h/2, y(i-1) + k2/2);
    k4 = h * f(x(i-1) + h, y(i-1) + k3);
    y(i) = y(i-1) + (k1 + 2*k2 + 2*k3 + k4)/6;
    x(i) = x(i-1) + h;
end

for i = 5:n+1
    x(i) = x(i-1) + h;
    f3 = f(x(i-1), y(i-1));
    f2 = f(x(i-2), y(i-2));
    f1 = f(x(i-3), y(i-3));
    f0 = f(x(i-4), y(i-4));
    ynew = y(i-1) + h*(55*f3 - 59*f2 + 37*f1 -9*f0) / 24;
    y(i) = y(i-1) + h*(9*f(x(i), ynew) + 19*f3 -5*f2 + f1)/24;
end

