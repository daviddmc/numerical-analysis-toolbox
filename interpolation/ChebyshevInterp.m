function yq = ChebyshevInterp( f, n ,a, b, xq)
% ChebyshevInterp    1D Chebyshev polynomial interpolation.
%   Yq = PolyInterp(F, N, A, B, Xq) interpolates to find Yq, the values 
%   of the underlying function F at the query points Xq using N-point
%   Chebyshev polynomial interpolation with in [A, B].
%
%   See also

%   Copyright 2017 Junshen Xu

x = (b+a)/2 + (b-a)/2 * cos((2*(1:n)-1) * (pi/2/n));
y = f(x);
yq = PolyInterp(x, y, xq);


end

