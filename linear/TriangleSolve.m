function X = TriangleSolve( A, B, shape)
%TriangleSolve   Solve a triangular linear system.
%   X = TriangleSolve(A, B) solves the linear system A*X = B, where A is a
%   lower triangle.
%
%   X = TriangleSolve(A, B, 'upper') solves the linear system A*X = B, 
%   where A is a upper triangle.
% 
%   See also

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(A);
CheckMultiplicationSize(A,[],B);
n = size(A, 1);
X = zeros(size(B));

if(~exist('shape','var'))
    shape = 'lower';
end

if(strcmpi(shape,'lower'))
    for ii = 1 : n
        X(ii, :) = (B(ii, :) - A(ii, 1:ii-1) * X(1:ii-1, :)) / A(ii,ii);
    end
elseif(strcmpi(shape, 'upper'))
    for ii = n : -1 : 1
        X(ii, :) = (B(ii, :) - A(ii, ii+1:end) * X(ii+1:end, :)) / A(ii,ii);
    end
else
    error('Shape flag must be ''upper'' or ''lower''.');
end

end

