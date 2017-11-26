function x = Bisection(f, a, b, maxIter, tolerance)
% bisection method for solving nonlinear equation f(x) = 0 in [a,b]
% input
% f : function
% a
% b
% maxIter : max number of iterations, default: 10
% tolerance : the tolerance of the method, default: 1e-6
% output
% x

if(~isscalar(a) || ~isreal(a) || ...
        ~isscalar(b) || ~isreal(b) || ...
        a >= b)
    error('a,b should be real numbers such that a < b');
end

if(~exist('maxIter','var') || isempty(maxIter))
    maxIter = 10;
end

if(~exist('tolerance','var') || isempty(tolerance))
    tolerance = 1e-6;
end

for k = 1 : maxIter
    m = (a+b)/2;
    fLeft = f(a);
    fMid = f(m);
    if(fLeft * fMid < 0)
        b = m;
    else
        a = m;
    end
    if(abs(b - a)<tolerance)
        break;
    end
end

x = (a + b) / 2;



