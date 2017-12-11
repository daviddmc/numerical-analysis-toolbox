function X = FullRankLS( A, B, method )
%LS Summary of this function goes here
%   Detailed explanation goes here

if(~ismatrix(A) || size(A,1) < size(A,2))
    error(' ');
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
    error('method');
end

