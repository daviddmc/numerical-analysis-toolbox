function [x, flag, iter, xs] = Brent(fun, interval, tol, maxIter)
% Brent   Find zero of single-variable function using Brent method.
%   X = Brent(FUN, INTERVAL) returns a zero of FUN within INTERVAL
%   using Brent method. FUN should be a function handle taking a scalar 
%   as input and INTERVAL should be a vector of length 2.
%
%   X = Brent(FUN, ITERVAL, TOL) specifies the tolerance of the method. 
%   If TOL is [] then Brent uses the default, 1e-6. 
%   
%   X = Brent(FUN, ITERVAL, TOL, MAXITER) specifies the maximum number 
%   of iterations. If MAXITER is [] then Brent uses the default,
%   max(10, ceil(log2(abs(INTERVAL(1)-INTERVAL(2)) / TOL)))
%
%   [X, FLAG] = Brent(FUN, ITERVAL, ...) also returns a convergence 
%   FLAG:
%    0 Brent converged to the desird tolerance TOL within MAXITER
%      iterations.
%    1 Brent iterated MAXITER times but did not converge.
%
%   [X, FLAG, ITER] = Brent(FUN, ITERVAL, ...) also returns the iteration 
%   number at which X was computed: 0 <= ITER <= MAXITER.
%
%   [X, FLAG, ITER, XS] = Brent(FUN, ITERVAL, ...) also returns a vector of
%   iteration result at each step, length(XS) = ITER, XS(end) = X;
%
%   See also

%   Copyright 2017 Junshen Xu

if(~isvector(interval) || length(interval) ~= 2)
    error('interval should be a vector of length 2.')
end

if(any(~isreal(interval)))
    error('interval should be real.')
end

a = interval(1);
b = interval(2);
fa = fun(a);
fb = fun(b);

if(fa * fb > 0)
    error('The function values at the interval endpoints must differ in sign.')
end

if(~exist('tol','var') || isempty(tol))
    tol = 1e-6;
end

if(~exist('maxIter','var') || isempty(maxIter))
    maxIter = max(10, ceil(log2((b-a)/tol)));
end

if(abs(fa) < abs(fb))
    temp = a; a = b; b = temp;
    temp = fa; fa = fb; fb = temp;
end

c = a;
fc = fun(c);
mflag = 1;
flag = 1;
if nargout > 3
    xs = zeros(maxIter, 1);
end
for iter = 1:maxIter
    if(fa ~= fc && fb ~= fc)
        s = a * fb * fc / (fa - fb) / (fa - fc);
        s = s + b * fa * fc / (fb - fa) / (fb - fc);
        s = s + c * fa * fb / (fc - fa) / (fc - fb);
    else
        s = b - fb * (b - a) / (fb - fa);
    end
    
    if((s - (3*a + b) / 4) * (s - b) > 0 || ...
            (mflag && abs(s-b) >= abs(b-c)/2) ||...
            (~mflag && abs(s-b) >= abs(c-d)/2) ||...
            (mflag && abs(b-c) < tol) ||...
            (~mflag && abs(c-d) < tol))
        s = (a+b)/2;
        mflag = 1;
    else
        mflag = 0;
    end
    fs = fun(s);
    d = c;
    c = b;
    if(fa*fs <= 0)
        b = s;
        fb = fun(b);
    else
        a = s;
        fa = fun(a);
    end
    
    if(abs(fa) < abs(fb))
        temp = a; a = b; b = temp;
        temp = fa; fa = fb; fb = temp;
    end
    
    if(nargout > 3)
        xs(iter) = s;
    end
    
    if(fs == 0 || abs(b - a) < tol)
        x = s;
        flag = 0;
        break
    end
    
end

if(nargout > 3)
    xs = xs(1:iter);
end


end

