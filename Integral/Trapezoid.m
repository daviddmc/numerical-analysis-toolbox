function I = Trapezoid( fun, a, b, n )
% Trapezoid     Numerical integration using trapezoid rule.
%   I = Trapezoid(FUN, A, B, N) approximates the integral of function FUN
%   from A to B using trapezod rule. The interval [A, B] will be divided
%   into N equispaced subintervals.
%
%   See also

%   Copyright 2017 Junshen Xu

h = (b - a) / n;
I = h * (fun(a) + fun(b) + 2 * sum(fun(a + h * (1:n-1)))) / 2;

end

