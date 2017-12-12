function [x, flag, iter] = Newtons(fun, jac, tol, maxIter, x0)
%NEWTONS Summary of this function goes here
%   Detailed explanation goes here
x = x0;
flag = 1;
for iter = 1:maxIter
    f = fun(x);
    J = jac(x);
    y = GaussianElimination(J, -f);
    x = x + y;
    if Norm(y) < tol
        flag = 0;
        break
    end
end

