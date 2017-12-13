function yInterp = HermiteRefine( tInterp, t0, y0, y1, f0, f1, h)
%HERMITEREFINE Summary of this function goes here
%   Detailed explanation goes here

s = (tInterp - t0)/h;
s2 = s .* s;
s3 = s .* s2;
slope = (y1 - y0) / h;
c = 3*slope - 2*f0 - f1;
d = f0 + f1 - 2*slope;
yInterp = y0(:,ones(size(tInterp))) + (h*d*s3 + h*c*s2 + h*f0*s);  
end

