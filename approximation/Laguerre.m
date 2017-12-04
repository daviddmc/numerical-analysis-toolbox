function L = Laguerre( n )
% Laguerre   Laguerre polynomials.
%   P = Laguerre(N) returns a vector P of length N+1 whose elements are the
%   coefficients of the Laguerre polynomial of order N in descending 
%   powers.
%
%   See also Chebyshev, Hermite, Legendre.

%   Copyright 2017 Junshen Xu

L0 = 1;
L1 = [-1 1];

if(n == 0)
    L = L0;
elseif(n == 1)
    L = L1;
else 
    for ii = 2 : n
        L = [0 (2*ii-1) * L1] - [L1 0];
        L(3:end) = L(3:end) - (ii-1)^2 * L0; 
        L0 = L1;
        L1 = L;
    end
end


