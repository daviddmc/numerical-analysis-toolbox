function [a, e] = NystromCoef( n )
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


g = [2, 0, 1/3, 1/3, 29/90, 14/45, 1139/3780, 41/140, 32377/113400];
m = length(g);

if n+1 <= m
    g = g(1:n+1);
else
    g = [g zeros(1, n+1-m)];
    r = 1 ./ (n+1:-1:1);
    for j = m+1:n+1
        g(j) = 1 - r(end-j+1:end-1)*g(1:j-1)';
    end
end

if n > 2
    e = g(end)/2;
else
    e = 1/6;
end

a = zeros(1,n);
c = [1 zeros(1,n)];
for j = 1 : n
    a(1:j) = a(1:j) + g(j) * c(1:j);
    c(2:(j+1)) = c(2:(j+1)) - c(1:j);
end

end