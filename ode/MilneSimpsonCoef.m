function [a, e] = MilneSimpsonCoef( n )
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


g = [2, -2];
m = length(g);

if n+1 <= m
    g = g(1:n+1);
else
    g = [g zeros(1, n+1-m)];
    r = 1 ./ (n+1:-1:1);
    for j = m+1:n+1
        g(j) = - r(end-j+1:end-1)*g(1:j-1)';
    end
end

% to do
if n > 3
    e = g(end)/2;
else
    e = -1/180;
end

a = zeros(1,n);
c = [1 zeros(1,n)];
for j = 1 : n
    a(1:j) = a(1:j) + g(j) * c(1:j);
    c(2:(j+1)) = c(2:(j+1)) - c(1:j);
end

end