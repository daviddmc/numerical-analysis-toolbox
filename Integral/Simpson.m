function I = Simpson( f, a, b, n )
%TRAPEZOID Summary of this function goes here
%   Detailed explanation goes here

h = (b - a) / n;
H = h * sum(f(a-h/2 + h * (1:n)));
T = Trapezoid(f, a, b, n);

I = T/3 + 2*H/3;
end

