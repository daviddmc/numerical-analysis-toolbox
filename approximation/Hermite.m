function H = Hermite( n)
% Hermite   Hermite polynomials.
%   P = Hermite(N) returns a vector P of length N+1 whose elements are the
%   coefficients of the Hermite polynomial of order N in descending powers.
%
%   See also Chebyshev, Legendre, Laguerre.

%   Copyright 2017 Junshen Xu

H0 = 1;
H1 = [2, 0];

if(n == 0)
    H = H0;
elseif(n == 1)
    H = H1;
else 
    for ii = 2 : n
        H = [2*H1 0];
        H(3:end) = H(3:end) - (ii-1)*2 * H0; 
        H0 = H1;
        H1 = H;
    end
end