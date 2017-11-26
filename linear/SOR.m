function X = SOR(A, B, omega ,X0, maxIter, tolerance)
% SOR

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

if(~exist('omega','var') || isempty(omega))
    omega = 1;
end

D = diag(diag(A));
invDL = Inverse(D + omega * tril(A, -1));
BSOR = invDL * ((1-omega)*D - omega * triu(A, 1)); 
FSOR = omega * invDL * B;

for k = 1 : maxIter
    X = BSOR * X0 + FSOR;
    if(max(abs(X(:) - X0(:))) < tolerance)
        break;
    end
    X0 = X;
end

