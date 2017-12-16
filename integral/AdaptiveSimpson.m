function I = AdaptiveSimpson( f, a, b, tol)
% AdaptiveSimpson    Numerical integration using adaptive Simpson method.
%   I = AdaptiveSimpson(FUN, A, B) approximates the integral of function 
%   FUN from A to B using adaptive Simpson method.
%
%   I = AdaptiveSimpson(FUN, A, B, TOL) specifies the tolerance. If TOL is
%   [], AdaptiveSimpson use the default, 1e-6.
%
%   See also

%   Copyright 2017 Junshen Xu

if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-6;
end

s = S(f, a, b);
m = (a+b)/2;
sLeft = S(f, a, m);
sRight = S(f, m, b);

if(abs(s - sLeft - sRight) < 10 * tol)
    I = sLeft + sRight;
else
    I = AdaptiveSimpson(f, a, m, tol/2) + ...
        AdaptiveSimpson(f, m, b, tol/2);
end

end

function s = S(f, a, b)
    s = (b - a) * (f(a) + 4 * f((a+b)/2) + f(b)) / 6;
end