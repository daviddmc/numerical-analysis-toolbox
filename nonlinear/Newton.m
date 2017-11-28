function [x, flag] = Newton( f, df, x0, maxIter, tolerance)
% Newton method

flag = 0;

if(~exist('x0','var') || isempty(x0))
    x0 = 0;
end

if(~exist('maxIter','var') || isempty(maxIter))
    maxIter = 10;
end

if(~exist('tolerance','var') || isempty(tolerance))
    tolerance = 1e-6;
end


for k = 1 : maxIter
    x = x0 - f(x0) / df(x0);
    if(abs(x - x0) < tolerance)
        flag = 1;
        break
    end
    x0 = x;
end