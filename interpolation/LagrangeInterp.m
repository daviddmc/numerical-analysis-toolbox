function yq = LagrangeInterp(X, y, Xq)
% 1-D polymonial interpolation using Lagrange method
% 

yq = zeros(size(Xq));
for i = 1 : length(yq)
    xq = Xq(i);
    for j = 1 : length(X)
        l = (xq - X) ./ (X(j) - X);
        yq(i) = yq(i) + prod(l([1:j-1 j+1:end])) * y(j);
    end
end