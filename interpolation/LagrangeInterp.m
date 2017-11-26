function L = LagrangeInterp(x, y)

    L = @(xq)lagrangeInterp(y, x, xq);

end

function Yq = lagrangeInterp(Y, X, Xq)
    Yq = zeros(size(Xq));
    for i = 1 : length(Yq)
        xq = Xq(i);
        for j = 1 : length(X)
            l = 1;
            for k = 1 : length(X)
                if(j ~= k)
                    l = l * (xq - X(k)) / (X(j) - X(k));
                end
            end
            Yq(i) = Yq(i) + l * Y(j);
        end
    end
end