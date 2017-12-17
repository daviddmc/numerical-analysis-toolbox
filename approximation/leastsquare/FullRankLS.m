function X = FullRankLS( A, B, method )
%FullRankLS    Full rank least square problem.
%   X = FullRankLS(A, B) solves the least square problem 
%       
%       argmin_x||A*x-B||_2
%
%   where A is a M-by-N matrix, M>=n and rank(A)=n.
%
%   X = FullRankLS(A, B, METHOD) specifies the method. The available
%   methods are:
%
%       QR     - (default) QR decomposition
%       normal - solving normal equation
%       SVD    - singular value decomposition
%
%   See also 

%   Copyright 2017 Junshen Xu

if(~ismatrix(A) || size(A,1) < size(A,2))
    error('A must be a matrix with full column rank.');
end

if ~exist('method', 'var')
    method = 'QR';
end

if(strcmp(method, 'normal'))
   X = SPDSolve(A' * A, A' * B);
elseif(strcmp(method, 'QR'))
   [Q, R] = QR(A, 'GramSchmidt');
   X = TriangleSolve(R, Q'*B, 'u');
elseif(strcmp(method, 'SVD'))
    [U, S, V] = SVD(A, 0);
    s = diag(S);
    X = V * (s .* (U' * B));
else
    error('method error');
end

