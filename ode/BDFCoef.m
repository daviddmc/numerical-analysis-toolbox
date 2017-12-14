function [a, e] = BDFCoef( n )
% AdamsCoef     Coefficient of Adams methods.
%   A = AdamsCoef(n) returns the coefficients in n-order explicit Adams
%   method:
%
%       Y(k+1) = Y(k) + h(A(1)*F(k) + A(2)*F(k-1) + ... + A(n)*F(k+1-n))
%   
%   [A, e] = AdamsCoef(n) also returns the coefficient of the main error in
%   n-order explicit Adams methods, so that
%
%       Y(x_k) - Y_k = h^(n+1)*e*Y^(k+1)(x_0)
%

if n > 6
    warning('For n > 6 the BDF-methods are unstable');
end

e = -1/(n+1);

a = zeros(1,n+1);
c = [1 zeros(1,n)];
for j = 1 : n
    c(2:(j+1)) = c(2:(j+1)) - c(1:j);
    a(1:j+1) = a(1:j+1) + c(1:j+1) / j;
end

end