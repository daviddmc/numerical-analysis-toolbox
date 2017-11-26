function X = JacobiIteration(A, B, X0, maxIter, tolerance)
% Jacobi iteration

CheckSquareMatrix(A);

if(~exist('maxIter', 'var') || isempty(maxIter))
    maxIter = 10;
end

if(~exist('X0','var') || isempty(X0))
    X0 = zeros(size(B));
end

CheckMultiplicationSize(A,X0, B);

if(~exist('tolerance', 'var') || isempty(tolerance))
    tolerance = 1e-6;
end

D = diag(A);
n = size(A, 1);
BJ = eye(n);
FJ = B;
for ii = 1 : n
    BJ(ii, :) = BJ(ii, :) - A(ii, :) / D(ii);
    FJ(ii, :) = B(ii, :) / D(ii);
end

for k = 1 : maxIter
    X = BJ * X0 + FJ;
    if(max(abs(X(:) - X0(:))) < tolerance)
        break;
    end
    X0 = X;
end

