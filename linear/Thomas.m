function X = Thomas(a, b, c, D)
% Thomas method for solving tridiagonal systems of equations
% [b1 c1
% [a1 b2 c2
% [   a2 b3 \\
% [      \\ \\ cn-1
% [        an-1  bn
% input
% A
% B
% output
% X

n = length(b);
if(length(a)~=n-1 || length(c)~=n-1)
    error(' ');
end

u = zeros(size(b));
l = zeros(size(a));
Y = zeros(size(D));
X = zeros(size(D));
u(1) = b(1);
for ii = 2:n
    l(ii-1) = a(ii - 1) / u(ii - 1);
    u(ii) = b(ii) - l(ii-1) * c(ii - 1);
end

Y(1,:) = D(1,:);
for ii = 2 : n
    Y(ii, :) = D(ii,:) - l(ii-1) * Y(ii-1,:);
end
X(n,:) = Y(n, :) / u(n);
for ii = n-1 : -1 : 1
    X(ii,:) = (Y(ii,:) - c(ii) * X(ii+1,:))/u(ii);
end






