function I = Romberg( fun, a, b, kmax, tol)
% Romberg     Numerical integration using Romber method.
%   I = Romberg(FUN, A, B) approximates the integral of function FUN
%   from A to B using Romberg method.
%
%   I = Romberg(FUN, A, B, Kmax) specifies the maximum of number of stages.
%   If Kmax is [], Romberg use the default, 5.
%
%   I = Romberg(FUN, A, B, Kmax, TOL) specifies the tolerance. If TOL is
%   [], Romberg use the default, 1e-6.
%
%   See also

%   Copyright 2017 Junshen Xu

if ~exist('kmax','var') || isempty(kmax)
    kmax = 5;
end

if ~exist('tol','var') || isempty(tol)
    tol = 1e-6;
end

T = zeros(kmax+1, 1);
T(1) = Trapezoid(fun, a, b, 1);
I = T(1);
for k = 1 : kmax
    T1 = Trapezoid(fun, a, b, 2^k);
    for jj = 1 : k+1
        T2 = (4^jj * T1 - T(jj)) / (4^jj - 1);
        T(jj) = T1;
        T1 = T2;
    end
    if(abs(I - T(k+1)) < tol)
        I = T(k+1);
        break
    end
    I = T(k+1);
end

end

