function x = Continuation( fun, jac, x0, N)
% Continuation   homotopy continuation method.
%   X = Continuation(FUN, JAC X0, N) starts at the column vector X0 and 
%   tries to solve the equations in FUN. FUN accepts input X and returns 
%   a column vector of equation values FUN evaluated at X. JAC accepts 
%   input X and returns the Jacobian matrix at X. N is the number of steps.
%
%   See also Newtons, Broyden.

%   Copyright 2017 Junshen Xu

h = 1/N;
b = -h*fun(x0);
x = x0;

for k = 1:N
    k1 = GaussianElimination(jac(x), b);
    k2 = GaussianElimination(jac(x + 0.5*k1), b);
    k3 = GaussianElimination(jac(x + 0.5*k2), b);
    k4 = GaussianElimination(jac(x + k3), b);
    x = x + (k1+2*k2+2*k3+k4) / 6;
end

end

