function [x, flag, iter, xs] = Bisection(fun, interval, tol, maxIter)
% Bisection   Find zero of single-variable function using bisection method.
%   X = Bisection(FUN, INTERVAL) returns a zero of FUN within INTERVAL
%   using bisection method. FUN should be a function handle taking a scalar 
%   as input and INTERVAL should be a vector of length 2.
%
%   X = Bisection(FUN, ITERVAL, TOL) specifies the tolerance of the method. 
%   If TOL is [] then Bisection uses the default, 1e-6. 
%   
%   X = Bisection(FUN, ITERVAL, TOL, MAXITER) specifies the maximum number 
%   of iterations. If MAXITER is [] then Steffensen uses the default,
%   max(10, ceil(log2(abs(INTERVAL(1)-INTERVAL(2)) / TOL)))
%
%   [X, FLAG] = Bisection(FUN, ITERVAL, ...) also returns a convergence 
%   FLAG:
%    0 Bisection converged to the desird tolerance TOL within MAXITER
%      iterations.
%    1 Bisection iterated MAXITER times but did not converge.
%
%   [X, FLAG, ITER] = Bisection(FUN, ITERVAL, ...) also returns the 
%   iteration number at which X was computed: 0 <= ITER <= MAXITER.
%
%   [X, FLAG, ITER, XS] = Bisection(FUN, ITERVAL, ...) also returns a 
%   vector of iteration result at each step, length(XS) = ITER, 
%   XS(end) = X;
%
%   See also

%   Copyright 2017 Junshen Xu

if(~isvector(interval) || length(interval) ~= 2)
    error('interval should be a vector of length 2.')
end

if(any(~isreal(interval)))
    error('interval should be real.')
end

a = min(interval);
b = max(interval);

if(fun(a) * fun(b) > 0)
    error('The function values at the interval endpoints must differ in sign.')
end

if(~exist('tol','var') || isempty(tol))
    tol = 1e-6;
end

if(~exist('maxIter','var') || isempty(maxIter))
    maxIter = max(10, ceil(log2((b-a)/tol)));
end

flag = 1;
if nargout > 3
    xs = zeros(maxIter, 1);
end
for iter = 1 : maxIter
    m = (a+b)/2;
    fLeft = fun(a);
    fMid = fun(m);
    if(nargout > 3)
        xs(iter) = m;
    end
    if(abs(fMid) < tol)
        flag = 0;
        break;
    end
    if(fLeft * fMid <= 0)
        b = m;
    else
        a = m;
    end
end

if(nargout > 3)
    xs = xs(1:iter);
end

x = m;



