function [ lambda, v ] = RayleighQuotientIteration( A, v0, maxIter, tolerance)
%RAYLEIGHQUOTIENTITERATION Summary of this function goes here
%   Detailed explanation goes here

CheckSquareMatrix(A);

if(~exist('maxIter', 'var') || isempty(maxIter))
    maxIter = 10;
end

if(~exist('v0','var') || isempty(v0))
    v0 = ones(size(A, 1), 1);
end

v = v0 / Norm(v0);

CheckMultiplicationSize(A, v0, []);

if(~exist('tolerance', 'var') || isempty(tolerance))
    tolerance = 1e-6;
end

n = size(A);
for k = 1 : maxIter
    lambda = v' * A * v;
    v = CholeskySolve(A - lambda * eye(n),v,'LDL');
    v = v / Norm(v);
    if(Norm(A*v - lambda * v) < tolerance)
        break
    end
end


end

