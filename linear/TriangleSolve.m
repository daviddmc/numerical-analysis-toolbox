function X = TriangleSolve( A, B, type)
%TRIANGLESOLVE Summary of this function goes here
%   Detailed explanation goes here

CheckSquareMatrix(A);
CheckMultiplicationSize(A,[],B);
n = size(A, 1);
X = zeros(size(B));

if(type == 'l')
    for ii = 1 : n
        X(ii, :) = (B(ii, :) - A(ii, 1:ii-1) * X(1:ii-1, :)) / A(ii,ii);
    end
elseif(type == 'u')
    for ii = n : -1 : 1
        X(ii, :) = (B(ii, :) - A(ii, ii+1:end) * X(ii+1:end, :)) / A(ii,ii);
    end
else
    error('type');
end

end

