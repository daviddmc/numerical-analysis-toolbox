function X = TriangleSolve( A, B, shape)
%TriangleSolve   Solve a triangular linear system.
%   X = TriangleSolve(A, B) solves the linear system A*X = B, where A is a
%   lower triangle.
%
%   X = TriangleSolve(A, B, 'upper') solves the linear system A*X = B, 
%   where A is a upper triangle.
% 
%   See also GaussElimination, SPDSolver, TridiagSolve.

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(A);
CheckMultiplicationSize(A,[],B);
n = size(A, 1);
X = zeros(size(B));

if(min(abs(diag(A))) < eps)
    warning('Matrix is singular to working precision.');
end

tol = 1e-14 * max(abs(diag(A)));

if(~exist('shape','var'))
    shape = 'lower';
end

if(strncmpi(shape,'lower', 1))
    for ii = 1 : n
        if(abs(A(ii, ii)) <= tol)
            X(ii, :) = 1;
        else
            X(ii, :) = (B(ii, :) - A(ii, 1:ii-1) * X(1:ii-1, :)) / A(ii,ii);
        end
    end
elseif(strncmpi(shape, 'upper', 1))
    for ii = n : -1 : 1
        if(abs(A(ii, ii)) <= tol)
            X(ii, :) = 1;
        else
            X(ii, :) = (B(ii, :) - A(ii, ii+1:end) * X(ii+1:end, :)) / A(ii,ii);
        end
    end
else
    error('Shape flag must be ''upper'' or ''lower''.');
end

end

