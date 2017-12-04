function DP = PolyDiff( P )
%PolyDiff   Differentiate polynomial.
%   PolyDiff(P) returns the derivative of the polynomial whose
%   coefficients are the elements of vector P.
%
%   See also 

%   Copyright 2017 Junshen Xu

if(~isvector(P))
    error('P should be a vector');
end

n = length(P) - 1;
DP = P(1:n) .* (n:-1:1);

end

