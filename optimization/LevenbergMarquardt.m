function [ x, flag, iter ] = LevenbergMarquardt( fun, jac, tol1, tol2, maxIter, x0, tau)
%LEVENBERGMARQUARDT Summary of this function goes here
%   Detailed explanation goes here

n = length(x0);
I = eye(n);
%tau = 1e-3;
x = x0;
v = 2;
J = jac(x);
f = fun(x);
F = f'*f / 2;
A = J'*J;
g = J'*f;
mu = tau * max(diag(A));
flag = 1;
if max(abs(g)) < tol1
    flag = 0;
    iter = 0;
    return
end
for iter = 1:maxIter
    h = SPDSolve(A+mu*I, -g);
    if Norm(h) < tol2*(Norm(x) + tol2)
        flag = 0;
        break
    end
    xNew = x + h;
    fNew = fun(xNew);
    FNew = fNew'*fNew / 2;
    rho = (F - FNew) / (h'*(mu*h-g)/2);
    if rho > 0
        x = xNew;
        f = fNew;
        F = FNew;
        J = jac(x);
        A = J'*J;
        g = J'*f;
        if max(abs(g)) < tol1
            flag = 0;
            break
        end
        mu = mu * max([1/3, 1-(2*rho-1)^3]);
        v = 2;
    else
        mu = mu * v;
        v = 2*v;
    end
        
end

end


