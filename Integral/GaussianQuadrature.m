function I = GaussianQuadrature( f, a, b, n ,method )
%GAUSSLEGENDRE Summary of this function goes here
%   Detailed explanation goes here


if(strcmp(method, 'Legendre'))
    p = Legendre(n + 1);
    pp = p(end:-1:1);
    pp = pp(1 : n+1) .* (n+1 : -1: 1);
    r = PolyRealRoot(p);
    A = 2./((1-r.^2) .* polyval(pp, r).^2);
    r = (a + b) / 2 + (b - a) / 2 * r;
    I = A * f(r)';
    I = (b - a) / 2 * I;
elseif(strcmp(method, 'Chebyshev'))
    r = 0 : n;
    r = cos(pi / (2*n+2) * (2*r + 1));
    r = (a + b) / 2 + (b - a) / 2 * r;
    I = sum(f(r));
    I = pi*(b - a) / (2*(n+1)) * I;
    
elseif(strcmp(method, 'Laguerre'))
    
elseif(strcmp(method, 'Hermite'))
    
else
    error('method');
end

