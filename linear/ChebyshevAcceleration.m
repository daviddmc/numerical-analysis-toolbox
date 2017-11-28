function X = ChebyshevAcceleration(R, C, rho,X0, maxIter, tolerance)
% ChebyshevAcceleration

CheckSquareMatrix(R);

if(~exist('maxIter', 'var') || isempty(maxIter))
    maxIter = 10;
end

if(~exist('X0','var') || isempty(X0))
    X0 = zeros(size(C));
end

CheckMultiplicationSize(R,X0, []);

if(~exist('tolerance', 'var') || isempty(tolerance))
    tolerance = 1e-6;
end

mu0 = 1;
mu1 = rho;
Y0 = X0;
Y1 = R * X0 + C;

for k = 2 : maxIter
    mu = 1 / (2 / (rho * mu1) - 1/ mu0);
    a = (2*mu / (rho * mu1));
    Y = R * Y1 * a - (mu / mu0) * Y0 + a * C;
    if(max(abs(Y(:) - Y1(:))) < tolerance)
        break;
    end
    Y0 = Y1;
    Y1 = Y;
    mu0 = mu1;
    mu1 = mu;
end

