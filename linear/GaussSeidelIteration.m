function X = GaussSeidelIteration(A, B, X0, maxIter, tolerance)
% Gauss - Seidel iteration

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

n = size(A, 1);
invDL = Inverse(tril(A, 0));
BG = eye(n) - invDL * A;
FG = invDL * B;

for k = 1 : maxIter
    X = BG * X0 + FG;
    if(max(abs(X(:) - X0(:))) < tolerance)
        break;
    end
    X0 = X;
end

