function I = AdaptiveSimpson( f, a, b, tol)
%ADAPTIVESIMPSON Summary of this function goes here
%   Detailed explanation goes here

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