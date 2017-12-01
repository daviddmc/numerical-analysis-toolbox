function yq = HermiteInterp( x, y, yp, xq, method)
%HERMITEINTERP Summary of this function goes here
%   Detailed explanation goes here

if(strcmpi(method, 'Lagrange'))
    yq = HermiteLagrange(x, y, yp, xq);
elseif(strcmpi(method, 'Newton'))
    yq = HermiteNewton(x, y, yp, xq);
end

end

function yq = HermiteLagrange(x, y, yp, Xq)

yq = zeros(size(Xq));
lp = 1 ./ bsxfun(@minus, x(:), x(:).');
for j = 1 : length(x)
    lp(j, j) = sum(lp(j, [1:j-1 j+1:end]));
end

for i = 1 : length(yq)
    xq = Xq(i);
    for j = 1 : length(x)
        l = (xq - x) ./ (x(j) - x);
        lsq = prod(l([1:j-1 j+1:end]))^2;
        yq(i) = yq(i) + ...
            lsq * (y(j) * (1 + 2 * lp(j, j) * (x(j) - xq)) + yp(j) * (xq - x(j)));
    end
end


end

function yq = HermiteNewton(x, y, yp, xq)

n = length(x);
nq = length(xq);

z = zeros(2*n,1);
z(1:2:end) = x;
z(2:2:end) = x;
Q = zeros(2*n);
Q(1:2:end, 1) = y;
Q(2:2:end, 1) = y;
Q(2:2:end, 2) = yp;

Q(3:2:end, 2) = (y(1:end-1) - y(2:end)) / (x(1:end-1) - x(2:end));

for ii = 3 : 2*n
    Q(ii : end, ii) = ...
        (Q(ii-1 : end - 1,ii-1) - Q( ii : end, ii-1)) ./...
        (z(1 : end - ii + 1) - z(ii : end));
end

Q = diag(Q);

A = bsxfun(@minus, xq(:) , [0 z(1:end-1).']);
A(:, 1) = 1;
A = cumprod(A, 2);

yq = A * Q;

end