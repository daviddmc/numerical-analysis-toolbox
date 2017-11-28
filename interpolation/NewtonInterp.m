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