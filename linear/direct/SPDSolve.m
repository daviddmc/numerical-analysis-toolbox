function X = SPDSolve(A, B)
% SPDSolve Solve semi-positive definite linear system
%   X = PDSolve(A, B) solves the positive definite linear system A*X=B 
%   using LDL dcomposition with pivoting. A is assumed to be 
%   symmetric semi-positive definite.
%
%   See also GaussElimination, TridiagSolve, TriangleSolve.

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(A);
CheckMultiplicationSize(A,[],B);
n = size(A, 1);

diagonalIdx = 1:n+1:n^2;
tol = max(Norm(A, inf) * 1e-14, 1e-14);
k = 1;
p = 1:n;
while k <= n
    [~, pivotK] = max(A(diagonalIdx(k:n)));
    pivotK = pivotK + k - 1;
    p([pivotK, k]) = p([k, pivotK]);
    A([pivotK, k],:) = A([k, pivotK],:);
    A(:, [pivotK, k]) = A(:,[k, pivotK]);
    A(k,k) = real(A(k,k));
    alpha = A(k, k);
    if(alpha < tol)
        break;
    end
    v = A(k+1:n, k);
    A(k+1:n, k) = v /alpha;
    A(k+1:n, k+1:n) = A(k+1:n, k+1:n) - v*v'/alpha;
    k = k+1;
end


X = B(p,:);
for i = 1:n
    X(i,:) = X(i,:) - A(i,1:i-1)*X(1:i-1,:);
end

for i = k-1:-1:1
    X(i,:) = X(i,:)/A(i, i) - A(i+1:n,i)'*X(i+1:n,:);
end

X(p,:) = X;

