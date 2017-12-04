function T = Chebyshev( n )
% Chebyshev   Chebyshev polynomials of the first kind.
%   P = Chebyshev(N) returns a vector P of length N+1 whose elements are 
%   the coefficients of the N-degree Chebyshev polynomial of the first kind
%   N in descending powers.
%
%   See also Legendre, Hermite, Laguerre.

%   Copyright 2017 Junshen Xu

T1 = 1;
T2 = [1, 0];

if(n == 0)
    T = T1;
elseif(n == 1)
    T = T2;
else 
    for ii = 2 : n
        T = [2*T2 0];
        T(3:end) = T(3:end) - T1; 
        T1 = T2;
        T2 = T;
    end
end

