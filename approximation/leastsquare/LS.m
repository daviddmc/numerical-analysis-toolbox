function X = LS( A, B, method)
%LS    least square problem.
%   X = LS(A, B) solves the least square problem 
%       
%       argmin_x||A*x-B||_2
%
%   where A is a M-by-N matrix, which may be rank-deficient. If A is
%   rank-deficient. X is the solution of the corresponding normal equaiton
%   that has the smallest norm.
%
%   X = LS(A, B, METHOD) specifies the method. The available
%   methods are:
%
%       QR     - (default) QR decomposition
%       SVD    - singular value decomposition
%
%   See also 

%   Copyright 2017 Junshen Xu

if ~exist('method','var')
    method = 'qr';
end

if strncmpi(method, 'svd', 1)
    [S,U,V] = SVD(A, 0);
    s = diag(S);
    tol = max(size(A)) * eps(max(abs(s)));
    r = sum(s > tol);
    V(:,r+1:end) = [];
    U(:,r+1:end) = [];
    s(r+1:end) = [];
    s = 1./s(:);
    X = V*(s.*(U'*B));
elseif strncmpi(method, 'qr', 1)
    [m,n] = size(A);
    k = size(B,2);
    [Q,R,perm] = QR(A,0);
    tol = abs(R(1)).*(max(m,n).^2).*eps;
    rankT0 = sum(abs(diag(R)) > tol);
    if rankT0 == n
        z = Q'*B;
        X(perm,1:k) = R \ z;
    else
        [P,S] = QR(R(1:rankT0,:)',0);
        z = Q(:,1:rankT0)'*B;
        X(perm,1:k) = P * (S' \ z);
    end
else
    error('method error')
end


end

