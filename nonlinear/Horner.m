function [y, yp] = Horner(p, x)
% Horner    Evaluate polynomial using Horner algorithm.
%   Y = Horner(P, X) returns the value of a polynomial P evaluated at X. P
%   is a vector of length N+1 whose elements are the coefficients of the
%   polynomial in descending powers.
%
%       Y = P(1)*X^N + P(2)*X^N(N-1) + ... + P(N)*X + P(N+1)
%   
%   If X is a matrix or vector, the polynomial is evaluated at all points
%   in X.
%
%   [Y, YP] = Horner(P, X) returns the value of a polynomial P and it's
%   gradient YP evaluated at X.
%
%       YP = N*P(1)*X^(N-1) + ... + 2*P(N-1)*X + P(N)
%
%   See also

%   Copyright 2017 Junshen Xu   


nc = length(p);
y = p(1);

if nargout == 2
    yp = p(1);
    for i = 2:nc-1
        y = x .* y + p(i);
        yp = x .* yp + y;
    end
    y = x .* y + p(end);
else
    for i=2:nc
        y = x .* y + p(i);
    end
end

end

