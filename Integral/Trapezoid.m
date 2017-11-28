function I = Trapezoid( f, a, b, n )
%TRAPEZOID Summary of this function goes here
%   Detailed explanation goes here

h = (b - a) / n;

I = h * (f(a) + f(b) + 2 * sum(f(a + h * (1:n-1)))) / 2;


end

