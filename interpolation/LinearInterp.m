function yq = LinearInterp(x, y, xq)
% LinearInterp    1D Linear interpolation.
%   Yq = LinearInterp(X, Y, Xq) interpolates to find Yq, the values of the 
%   underlying function F(X) at the query points Xq using linear
%   interpolation.
%
%   See also

%   Copyright 2017 Junshen Xu

    [x, idx] = sort(x(:));
    y = y(idx);
    
    yq = zeros(size(xq));
    i = 1;
    while(i <= length(yq) && xq(i) < x(1))
        yq(i) = y(1);
        i = i + 1;
    end
    j = 2;
    while(j <= length(y) )
        while(i <= length(yq) && xq(i) < x(j))
            yq(i) = y(j-1)*(xq(i) - x(j))/(x(j-1) - x(j)) + y(j)*(xq(i) - x(j-1))/(x(j)-x(j-1));
            i = i + 1;
        end
        j = j + 1;
    end
    
    while(i <= length(yq))
        yq(i) = y(end);
        i = i + 1;
    end

end

