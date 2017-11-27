function [m, u] = InversePower( A, v0, shift, maxIter, tolerance )
% power method to calculate the largest eigenvalue

CheckSquareMatrix(A);

if(~exist('maxIter', 'var') || isempty(maxIter))
    maxIter = 10;
end

if(~exist('v0','var') || isempty(v0))
    v0 = ones(size(A, 1), 1);
end

CheckMultiplicationSize(A, v0, []);

if(~exist('tolerance', 'var') || isempty(tolerance))
    tolerance = 1e-6;
end

if(~exist('shift', 'var') || isempty(shift))
    shift = 0;
end

u = v0;
[~, idx] = max(abs(u));
m0 = u(idx);
n = size(A, 1);
for k = 1 : maxIter
    v = GaussianElimination(A - shift * eye(n), u);
    [~, idx] = max(abs(v));
    m = v(idx);
    u = v / m;
    if(abs(m - m0) < tolerance)
        break
    end
    m0 = m;
end

