function [x, flag] = Newton( f, df, x0, maxIter1, maxIter2, tolerance)
% Newton method

flag = 0;

if(~exist('x0','var') || isempty(x0))
    x0 = 0;
end

if(~exist('maxIter1','var') || isempty(maxIter1))
    maxIter1 = 10;
end

if(~exist('maxIter2','var') || isempty(maxIter2))
    maxIter2 = 10;
end

if(~exist('tolerance','var') || isempty(tolerance))
    tolerance = 1e-6;
end

if(df(x0) == 0)
    x0 = x0 + 10 * tolerance;
end

for k = 1 : maxIter1
    fx = f(x0);
    dfx = df(x0);
    delta = fx / dfx;
    x = x0 - delta;
    if(abs(delta) < tolerance || (abs(fx) < tolerance && abs(dfx) < tolerance))
        flag = 1;
        break
    end
    x0 = x;
end

if ~flag
    for k = 1 : maxIter2
        fx = f(x0);
        dfx = df(x0);
        delta = fx * log(abs(fx)) / dfx / (log(abs(fx)) - log(abs(dfx)));
        if(isnan(delta))
            flag = 1;
        else
            x = x0 - delta;
        end
        if(abs(delta) < tolerance || (abs(fx) < tolerance && abs(dfx) < tolerance))
            flag = 1;
        end
        if(flag)
            break;
        end
        x0 = x;
    end
end