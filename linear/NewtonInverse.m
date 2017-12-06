function [X, flag, iter] = NewtonInverse(A, tol, maxIter)
% NewtonInverse     pseudo inverse estimated by Newton iteration.
%   X = NewtonInverse(A) attempts to fine the pseudo inverse of matrix A 
%   using Newton method. If A is of size M-by-N, then X is N-by-M; 
%
%   X = NewtonInverse(A, TOL) specifies the tolerance of the method. If TOL 
%   is [] then NewtonInverse uses the default, 1e-6.
%
%   X = NewtonInverse(A, TOL, MAXITER) specifies the maximum number of
%   iterations. If MAXITER is [] then NewtonInverse uses the default, 
%   max(N, 20).
%
%   [X, FLAG] = NewtonInverse(A, ...) also returns a convergence FLAG:
%    0 NewtonInverse converged to the desird tolerance TOL within MAXITER
%      iterations.
%    1 NewtonInverse iterated MAXITER times but did not converge.
%    2 Input matrix is zero.
%
%   [X, FLAG, ITER] = NewtonInverse(A, ...) also returns the iteration
%   number at which X was computed: 0 <= ITER <= MAXITER.
%
%   See also

%   Reference
%   [1] S? T, derstr?, Stewart G W. On the Numerical Properties of an 
%   Iterative Method for Computing the Moore¨CPenrose Generalized Inverse. 
%   Siam Journal on Numerical Analysis, 1974, 11(1):61-74.
%   [2] Pan, Victor, and R. Schreiber. An improved Newton interaction for 
%   the generalized inverse of a Matrix, with applications. Siam, 1991.
%
%   Copyright 2017 Junshen Xu

if(~ismatrix(A))
    error('A should be a matrix');
end

n = size(A, 2);
if(~exist('maxIter', 'var') || isempty(maxIter))
    maxIter = max(20,n);
end

if(~exist('tol', 'var') || isempty(tol))
    tol = 1e-6;
end

alpha0 = 1 / (Norm(A, 1) * Norm(A, inf));
if(alpha0 == 0 || isnan(alpha0) || isinf(alpha0))
    flag = 2;
    iter = 0;
    return
end
tol = tol * alpha0;
X = alpha0 * A';
T = X*A;
T2valid = 0;
twoI = 2*eye(size(T));
i = 1;
flag = 1;
for iter = 1 : maxIter
    % check criteria
    if(iter > 1)
        if(trace(T) >= i)
            if(Norm(Xnew - X, 'fro') < tol * pow2(iter))
                flag = 0;
                break;
            else
                if i < n
                    i = i + 1;
                end
            end
        end
        X = Xnew;
    end
    % Newton step
    Xnew = (twoI - T) * X;
    if T2valid
        T = 2*T - T2;
    else
        T = Xnew*A;
    end
    T2valid = 0;
    if trace(T) >= n-0.5
        continue;
    end
    % test for small change
    T2 = T*T;
    delta = Norm(T-T2,'fro');
    if delta >= 0.25
        T2valid = 1;
        continue;
    end
    % use acceleration
    rho = 0.5 - sqrt(0.25 - delta);
    Xnew = 1/rho * (T2 - (2+rho) * T + (0.5+rho)*twoI) * Xnew;
    T = Xnew*A;
end

end

