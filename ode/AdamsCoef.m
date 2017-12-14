function [a, e] = AdamsCoef( n , implicit)
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

if exist('implicit','var') && strcmpi(implicit, 'implicit')
    a = 0;
    g = sym([1, -0.5, -1/12, -1/24, -19/720, -3/160, -863/60480, -275/24192, -33953/3628800]);
else
    g = [1, 0.5, 5/12, 3/8, 251/720, 95/288, 19087/60480, 5257/17280, 1070017/3628800];
    a = 1;
end
m = length(g);

if n+1 <= m
    g = g(1:n+1);
else
    g = [g zeros(1, n+1-m)];
    r = 1 ./ (n+1:-1:1);
    for j = m+1:n+1
        g(j) = a - r(end-j+1:end-1)*g(1:j-1)';
    end
end

e = g(end);

a = zeros(1,n);
c = [1 zeros(1,n)];
for j = 1 : n
    a(1:j) = a(1:j) + g(j) * c(1:j);
    c(2:(j+1)) = c(2:(j+1)) - c(1:j);
end

end