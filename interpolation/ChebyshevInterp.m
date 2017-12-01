function yq = ChebyshevInterp( f, n ,a, b, xq)
%CHEBYSHEVINTERP Summary of this function goes here
%   Detailed explanation goes here

x = (b+a)/2 + (b-a)/2 * cos((2*(1:n)-1) * (pi/2/n));
y = f(x);
yq = PolyInterp(x, y, xq);


end

