function yq = PolyInterp(x, y, xq, method)
%POLYINTERP Summary of this function goes here
%   Detailed explanation goes here

if(~exist('method', 'var') || isempty(method))
    method = 'newton';
end

if(strcmpi(method, 'newton'))
    yq = NewtonInterp(x, y, xq);
elseif(strcmpi(method, 'lagrange'))
    yq = LagrangeInterp(x, y, xq);
elseif(strcmpi(method, 'neville'))
    yq = Neville(x, y, xq);
else
    error(' ');
end

end

function yq = LagrangeInterp(X, y, Xq)
% 1-D polymonial interpolation using Lagrange method
% 

yq = zeros(size(Xq));
for i = 1 : length(yq)
    xq = Xq(i);
    for j = 1 : length(X)
        l = (xq - X) ./ (X(j) - X);
        yq(i) = yq(i) + prod(l([1:j-1 j+1:end])) * y(j);
    end
end
end

function yq = Neville( x, y, Xq)
%NEVILLE Summary of this function goes here
%   Detailed explanation goes here

%epsilon = 1e-6;
n = length(x);
Q = zeros(n);
Q(:, 1) = y;
yq = zeros(size(Xq));
for k = 1 : length(yq)
    xq = Xq(k);
    for i = 1 : n - 1
        for j = 1 : i
            Q(i+1, j+1) = ((xq - x(i-j+1)) * Q(i+1, j) - (xq - x(i+1)) * Q(i, j)) / (x(i+1) - x(i-j+1));
        end
    end
    yq(k) = Q(n ,n);
end
end

function yq = NewtonInterp(x, y, xq)
% 1-D polymonial interpolation using Newton method
% 

n = length(x);
D = zeros(n);
D(1,:) = y;

for ii = 2 : n
    D(ii, 1 : end - ii + 1) = ...
        (D(ii-1, 1 : end - ii + 1) - D(ii-1, 2 : end - ii + 2)) ./...
        (x(1 : end - ii + 1) - x(ii : end));
end
D = D(:,1);

if(iscolumn(x))
    x = x';
end

A = bsxfun(@minus, xq(:) , [0 x(1:end-1)]);
A(:, 1) = 1;
A = cumprod(A, 2);

yq = A * D;
end

