function [x, flag, iter, xs] = Secant(fun, x0, x1, tol, maxIter)
% Secant   Find zero of single-variable function using secant method.
%   X = Secant(FUN, X0, X1) returns a zero of FUN within X0, X1 as the
%   first two points using Secant method. FUN should be a function handle 
%   taking a scalar as input and X0, X1 should be scalar.
%
%   X = Secant(FUN, X0, X1, TOL) specifies the tolerance of the method. 
%   If TOL is [] then Secant uses the default, 1e-6. 
%   
%   X = Secant(FUN, X0, X1, TOL, MAXITER) specifies the maximum number 
%   of iterations. If MAXITER is [] then Steffensen uses the default, 10
%
%   [X, FLAG] = Secant(FUN, X0, X1, ...) also returns a convergence 
%   FLAG:
%    0 Secant converged to the desird tolerance TOL within MAXITER
%      iterations.
%    1 Secant iterated MAXITER times but did not converge.
%
%   [X, FLAG, ITER] = Secant(FUN, X0, X1, ...) also returns the iteration
%   number at which X was computed: 0 <= ITER <= MAXITER.
%
%   [X, FLAG, ITER, XS] = Secant(FUN, X0, X1, ...) also returns a vector of
%   iteration result at each step, length(XS) = ITER, XS(end) = X;
%
%   See also

%   Copyright 2017 Junshen Xu

if(~exist('maxIter','var') || isempty(maxIter))
    maxIter = 10;
end

if(~exist('tol','var') || isempty(tol))
    tol = 1e-6;
end

flag = 1;
if nargout > 3
    xs = zeros(maxIter, 1);
end
for iter = 1 : maxIter
    y1 = fun(x1);
    y0 = fun(x0);
    x = x1 - y1 * (x1 - x0) / (y1 - y0);
    if(abs(x - x1) < tol)
        flag = 0;
        break
    end
    x0 = x1;
    x1 = x;
    if(nargout > 3)
        xs(iter) = x;
    end
end