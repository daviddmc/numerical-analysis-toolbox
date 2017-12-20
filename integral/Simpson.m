function I = Simpson( f, a, b, n )
% Simpson     Numerical integration using Simpson rule.
%   I = Simpson(FUN, A, B, N) approximates the integral of function FUN
%   from A to B using Simpson rule. The interval [A, B] will be divided
%   into N equispaced subintervals.
%
%   See also

%   Copyright 2017 Junshen Xu

h = (b - a) / n;
H = h * sum(f(a-h/2 + h * (1:n)));
T = Trapezoid(f, a, b, n);

I = T/3 + 2*H/3;
end

