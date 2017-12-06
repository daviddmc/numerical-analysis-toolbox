function [M, U, flag] = EigenDeflation( A, tol, maxIter, symString)
% EigenDeflation     Find eigenvalues using iteration and deflation method.
%   M = EigenDeflation(A) tries to find all eigenvalues of matrix A using
%   power iteration and Wielandt deflation.
%
%   M = EigenDeflation(A, TOL) specifies the tolerance of the method. If 
%   TOL is [] then EigenDeflation uses the default, 1e-6.
%
%   M = EigenDeflation(A, TOL, MAXITER) specifies the maximum number of
%   iterations in each subprocess. If MAXITER is [] then EigenDeflation 
%   uses the default, max(N, 20).
%
%   M = EigenDeflation(A, TOL, MAXITER, 'symmetric') assumes that A is
%   symmetric(Hermitian) matrix. symmetric power iteration and Hotelling
%   deflation will be used.
%
%   [M, U] = EigenDeflation(A, ...) also returns the corresponding 
%   eigenvectors.
%
%   [M, U, FLAG] = EigenDeflation(A, ...) also returns a convergence FLAG:
%    0 EigenDeflation found all eigenvalue of A to the desird tolerance TOL 
%      within MAXITER iterations.
%    1 EigenDeflation is stopped for one of the subproccess iterated 
%      MAXITER times but did not converge.
%
%   See also PowerIter.

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(A);
n = size(A, 1);
if(~exist('maxIter', 'var') || isempty(maxIter))
    maxIter = max(20, n);
end

if(~exist('tol', 'var') || isempty(tol))
    tol = 1e-6;
end

M = zeros(n, 1);
flag = 0;
B = A;

if(exist('symString','var') && strncmpi(symString,'sym',1))
    U = zeros(n);
    for k = 1:n
        [M(k), u, subflag] = PowerIter(B, tol, maxIter, [], 's');
        U(:, k) = u;
        if subflag
            [M(k), u, subflag] = RayleighQuotientIter(A, tol, maxIter, u);
            if subflag
                flag = 1;
                break;
            end
        end
        B = B - M(k)/(u'*u)*(u*u');
    end
    M = M(1:k);
    U = U(:,1:k);
else
    for k = 1:n
        [M(k), u, subflag] = PowerIter(B, tol, maxIter);
        if subflag
            flag = 1;
            break;
        else
            [m, idx] = max(u);
            B = B - u/m * B(idx,:);
            B = B([1:idx-1 idx+1:end], [1:idx-1 idx+1:end]);
        end
    end
    M = M(1:k);

    if nargout > 1
        U = zeros(n, k);
        for k = 1:k
            [M(k), U(:,k)] = InverseIter(A, M(k), tol, maxIter);
        end
    else
        for k = 1:k
            M(k) = InverseIter(A, M(k), tol, maxIter);
        end
    end
    
end






