function yq = NearestInterp(x, y, xq)
    
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
            if(xq(i) - x(j-1) > x(j) - xq(i))
                yq(i) = y(j);
            else
                yq(i) = y(j-1);
            end
            %yq(i) = y(j-1)*(xq(i) - x(j))/(x(j-1) - x(j)) + y(j)*(xq(i) - x(j-1))/(x(j)-x(j-1));
            i = i + 1;
        end
        j = j + 1;
    end
    
    while(i <= length(yq))
        yq(i) = y(end);
        i = i + 1;
    end

end

