function [x, flag, iter, res] = GMRES( A, b ,tol, MaxIter, x0)
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

x = x0;
delta = max(tol, tol * Norm(b));
flag = 1;

for iter = 1 : MaxIter
    if 0
        r = M \ (b - A*x);
    else
        r = b - A*x;
    end
    rNorm = sqrt(r'*r);
    if flagRecord
        res(iter) = rNorm;
    end
    if rNorm < delta
        flag = 0;
        break
    end
    v = r / rNorm;
    s = rNorm*e1;
    for i = 1:m
        if 0
            w = M\(A*vi);
        else
            w = A*vi;
        end
        for k = 1:i
            h(k, i) = w'*vk;
            w = w - h(k,i)*vk;
        end
        h(i+1, i) = Norm(w);
        vi+1 = w / h(i+1, i);
        apply J
        get J
        s = J*s;
    
    
    
end

if flagRecord
    res = res(1:iter);
end

end

