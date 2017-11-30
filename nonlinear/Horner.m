function [y, yp] = Horner( p, x)
%HORNER Summary of this function goes here
%   Detailed explanation goes here


nc = length(p);
y = p(1);

if nargout == 2
    yp = p(1);
    for i = 2:nc-1
        y = x .* y + p(i);
        yp = x .* yp + y;
    end
    y = x .* y + p(end);
else
    for i=2:nc
        y = x .* y + p(i);
    end
end

end

