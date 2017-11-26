function X = Thomas(a, b, c, D)
% Thomas method for solving tridiagonal systems of equations
% [b1 c1
% [a2 b2 c2
% [   a3 b3 \\
% [      \\ \\ cn-1
% [         an  bn
% input
% A
% B
% output
% X

u = zeros(size(a));
l = zeros(size(a));
Y = zeros(size(D));
X = zeros(size(D));
u(1) = b(1);
for ii = 2:n
    l(ii) = a(ii) / u(ii - 1);
    u(ii) = b(ii) - l(ii) * c(ii - 1);
end

Y(1,:) = D(1,:);
for ii = 2 : n
    Y(ii, :) = D(ii,:) - l(ii) * Y(ii-1,:);
end
X(n,:) = Y(n, :) / u(n);
for ii = n-1 : -1 : 1
    X(ii,:) = (Y(ii,:) - c(ii) * X(ii+1,:))/u(ii);
end






