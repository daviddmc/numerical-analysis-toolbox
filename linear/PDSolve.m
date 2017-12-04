function X = PDSolve(A, B, method)
% PDSolve Solve positive definite linear system
%   X = PDSolve(A, B) solves the positive definite linear system A*X=B 
%   using LDL dcomposition. Since A is assumed to be symmetric positive 
%   definite, only the diagnoal and lower triangle of A will be used. 

CheckSquareMatrix(A);
CheckMultiplicationSize(A,[],B);
n = size(A, 1);
Y = zeros(size(B));
X = zeros(size(B));

if(~exist('method','var') || isempty(method))
    method = 'LDL';
end

if(strcmp(method, 'LL'))
    L = Cholesky(A);
    for ii = 1:n
        Y(ii,:) = (B(ii,:) - L(ii,1:ii-1) * Y(1:ii-1,:)) / L(ii,ii);
    end
    for ii = n : -1 : 1
        X(ii, :) = (Y(ii, :) - L(ii+1:n,ii)' * X(ii+1:n,:)) / L(ii, ii);
    end
elseif(strcmp(method, 'LDL'))
    [L, D] = Cholesky(A);
    for ii = 1:n
        Y(ii, :) = B(ii, :) - L(ii, 1:ii-1) * Y(1:ii-1,:);
    end
    for ii = n : -1 : 1
        X(ii, :) = Y(ii,:)/D(ii) - L(ii+1:n,ii)' * X(ii+1:n,:);
    end
else
    error('method');
end


end

