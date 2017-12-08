function P = Legendre( n )
% Legendre   Legendre polynomials.
%   P = Legendre(N) returns a vector P of length N+1 whose elements are the
%   coefficients of the Legendre polynomial of order N in descending
%   powers.
%
%   See also Chebyshev, Hermite, Laguerre.

%   Copyright 2017 Junshen Xu
 
P0 = 1;
P1 = [1, 0];

if(n == 0)
    P = P0;
elseif(n==1)
    P = P1;
else
    for ii = 2 : n
        P = [(2*ii - 1)*P1 0];
        P(3 : end) = P(3 : end) - (ii-1) * P0;
        P = P / ii;
        P0 = P1;
        P1 = P;
    end
       
end

