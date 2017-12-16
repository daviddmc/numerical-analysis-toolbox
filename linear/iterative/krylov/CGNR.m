function [x, flag, iter, res] = CGNR( A, b ,tol, MaxIter, x0)
% conjugate gradient method
% solving A * x = b
% input
% A
% b
% MaxIter : max number of iteration, default: size of A
% x0 : initial value of x
% tolerance : specifies the tolerance of the method, default: 1e-6
% output
% X

CheckSquareMatrix(A);
if(~isvector(b))
    error('b should be a vector');
end
if(~iscolumn(b))
    b = b(:);
end


if(~exist('MaxIter', 'var') || isempty(MaxIter))
    MaxIter = max(size(A,1), 20);
end

if(~exist('x0','var') || isempty(x0))
    x0 = zeros(size(A,1), 1);
end
if(~isvector(x0))
    error('x0 should be a vector');
end
if(~iscolumn(x0))
    x0 = x0(:);
end
CheckMultiplicationSize(A,x0, b);

if(~exist('tol', 'var') || isempty(tol))
    tol = 1e-6;
end

flagRecord = nargout > 3;
if flagRecord
    res = zeros(maxIter, 1);
end

delta = max(tol, tol * Norm(b));
delta2 = delta * delta;
flag = 1;

x = x0;
r = b - A * x;
z = A'*r;
rhoc = z'*z;
for iter = 1 : MaxIter
    if flagRecord
        res(iter) = sqrt(rhoc);
    end
    if rhoc < delta2
        flag = 0;
        break
    end
    if iter == 1
        p = z;
    else
        tau = rhoc / rho;
        p = z + tau * p;
    end
    w = A*p;
    mu = rhoc/(w'*w);
    x = x + mu*p;
    r = r - mu*w;
    z = A'*r;
    rho = rhoc;
    rhoc = z'*z;
end

if flagRecord
    res = res(1:iter);
end

end

