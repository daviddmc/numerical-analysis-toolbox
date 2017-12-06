function [m, u, flag, iter] = InverseIter( A, shift, tol,maxIter, u0)
% InverseIter     Inverse Power iteration method.
%   M = InverseIter(A) tries to find the smallest (in absolute value) 
%   eigenvalue of matrix A.
%
%   M = InverseIter(A, SHIFT) tries to find the eigenvalue of A that is
%   closest to SHIFT (in absolute value) If SHIFT is [] then InverseIter 
%   uses the default, 0.
%
%   M = InverseIter(A, SHIFT, TOL) specifies the tolerance of the method. 
%   If TOL is [] then PowerIter uses the default, 1e-6.
%
%   M = InverseIter(A, SHIFT, TOL, MAXITER) specifies the maximum number of
%   iterations. If MAXITER is [] then InverseIter uses the default, 
%   max(N, 20).
%
%   M = PowerIter(A, SHIFT, TOL, MAXITER, U0) specifies the initial guess 
%   of the eigenvector of the dominant eigenvalue. If U0 is [] then 
%   InverseIter uses the default, an all one vector. If U0 is specified
%   while SHIFT is [], SHIFT will be set to U0'*A*U0/(U0'*U0). 
%
%   [M, U] = PowerIter(A, SHIFT ...) also returns the corresponding 
%   eigenvector.
%
%   [M, U, FLAG] = PowerIter(A, SHIFT, ...) also returns a convergence 
%   FLAG:
%    0 InverseIter converged to the desird tolerance TOL within MAXITER
%      iterations.
%    1 InverseIter iterated MAXITER times but did not converge.
%    2 SHIFT is a eigenvalue of A.
%
%   [M, U, FLAG, ITER] = PowerIter(A, SHIFT ...) also returns the iteration
%   number at which M was computed: 0 <= ITER <= MAXITER.
%
%   See also PowerIter.

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(A);
n = size(A, 1);
if(~exist('maxIter', 'var') || isempty(maxIter))
    maxIter = max(20, n);
end

if(exist('u0','var') && isempty(shift))
    shift = (u0'*A*u0) / (u0'*u0);
end

if(~exist('u0','var') || isempty(u0))
    u0 = ones(size(A, 1), 1);
end

CheckMultiplicationSize(A, u0, []);

if(~exist('tol', 'var') || isempty(tol))
    tol = 1e-6;
end

if(~exist('shift', 'var') || isempty(shift))
    shift = 0;
end

[L, U, P] = LU(A - shift * eye(n));
if(min(abs(diag(U))) < 1e-14*max(abs(diag(U))))
    flag = 2;
    m = shift;
    iter = 0;
    u = 0*u0;
    u = TriangleSolve(L, P*u, 'lower');
    u = TriangleSolve(U, u, 'upper');
    return 
end
permuteIdx = zeros(1,n);
for k = 1 : n
    permuteIdx(k) = find(P(:,k));
end

if(all(permuteIdx == 1:n))
    permuteFlag = 0;
else
    permuteFlag = 1;
end

u = u0;
v = u0;
flag = 1;
for iter = 1 : maxIter
    if(permuteFlag)
        v(permuteIdx) = u;
    end
    v = TriangleSolve(L, v, 'lower');
    v = TriangleSolve(U, v, 'upper');  
    [~, idx] = max(abs(v));
    m = v(idx);
    if(max(abs(v - m*u)) < tol)
        flag = 0;
        break
    end
    u = v / m;
end

m = 1/m + shift;

