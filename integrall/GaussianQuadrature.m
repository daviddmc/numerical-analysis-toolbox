function I = GaussianQuadrature( f, a, b, n ,method )
%GAUSSLEGENDRE Summary of this function goes here
%   Detailed explanation goes here

if(strcmp(method, 'Legendre'))
    p = Legendre(n + 1);
    dp = PolyDiff(p);
    r = PolyRealRoot(p);
    A = 2./((1-r.^2) .* polyval(dp, r).^2);
    r = (a + b) / 2 + (b - a) / 2 * r;
    I = A * f(r)';
    I = (b - a) / 2 * I;
elseif(strcmp(method, 'Chebyshev'))
    r = 0 : n;
    r = cos(pi / (2*n+2) * (2*r + 1));
    I = sum(f(r));
    I = pi / (n+1) * I;
elseif(strcmp(method, 'Laguerre'))
    p = Laguerre(n + 1);
    dp = PolyDiff(p);
    r = PolyRealRoot(p);
    A = factorial(n+1)^2 ./(r .* polyval(dp, r).^2);
    I = A * f(r)';
elseif(strcmp(method, 'Hermite'))
    p = Hermite(n + 1);
    dp = PolyDiff(p);
    r = PolyRealRoot(p);
    A = factorial(n+1) * 2^(n+2) * sqrt(pi) ./ polyval(dp, r).^2;
    I = A * f(r)';
else
    error('method');
end

