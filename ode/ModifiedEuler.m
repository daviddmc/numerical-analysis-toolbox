function [x, y] = ModifiedEuler( f, a, b, y0, n)
%EULER Summary of this function goes here
%   Detailed explanation goes here

h = (b - a) / n;
x = zeros(n+1, 1);
y = zeros(n+1, 1);
x(1) = a;
y(1) = y0;

for ii = 1:n
    x(ii+1) = x(ii) + h;
    y1 = y(ii) + h * f(x(ii), y(ii));
    y(ii+1) = y(ii) + h/2 * (f(x(ii), y(ii)) + f(x(ii+1), y1));
end

