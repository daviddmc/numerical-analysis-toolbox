function yq = CubicSplineInterp( x, y, xq, boundary)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = length(x);
if(~exist('boundary', 'var') || isempty(boundary))
    if(n < 4)
        boundary = 'natural';
    else
        boundary = 'notaknot';
    end
end

if(n < 2)
    error('Too less points.');
elseif(n < 4 && strcmpi(boundary,'notaknot'))
    warning('Too less points for not-a-knot end conditions, use natural end condition instead.');
    boundary = 'natural';
end

[x, idx] = sort(x);
y = y(idx);

x = x(:);
y = y(:);
h = x(2:end) - x(1:end-1);

t = (y(2:end) - y(1:end-1)) ./ h;
alpha = t(2:end) - t(1:end-1);
alpha = 3 * alpha;

if(strcmpi(boundary, 'notaknot'))
    c = TridiagSolve([h(1:end-1); -3*(h(end)+h(end-1))*(1+2*h(end)/h(end-1))],...
        [3*(h(2) - h(1)/h(2)*h(1)); 2*(h(1:end-1)+h(2:end)); 3*(h(end-1) - h(end)/h(end-1)*h(end))],...
        [-3*(h(1)+h(2))*(1+2*h(1)/h(2));h(2:end)],...
        [-3*h(1)/h(2)*alpha(1); alpha ;-3*h(end)/h(end-1)*alpha(end)]);
elseif(strcmpi(boundary, 'natural'))
    c = TirdiagSolve([h(1:end-1);0], [1; 2*(h(1:end-1)+h(2:end)); 1] ,[0;h(2:end)], [0; alpha ;0]);
end

b = t - h .* (c(2:end) + 2*c(1:end-1)) / 3;
d = (c(2:end) - c(1:end-1)) ./ (3 * h);

[xq, idx] = sort(xq);
yq = zeros(size(xq));
m = length(xq);
i = 1;
j = 1;

while(xq(i) <= x(1))
    yq(i) = y(1) + b(1)*(xq(i) - x(1)) + c(1)*(xq(i) - x(1))^2 + d(1)*(xq(i)-x(1))^3;
    i = i + 1;
end

while(j < n && i <= m)
   while(i <= m && xq(i) <= x(j + 1))
       yq(i) = y(j) + b(j)*(xq(i) - x(j)) + c(j)*(xq(i) - x(j))^2 + d(j)*(xq(i)-x(j))^3;
       i = i + 1;
   end
   j = j + 1;
end

while(i<=m)
    yq(i) = y(n-1) + b(end)*(xq(i) - x(n-1)) + c(end-1)*(xq(i) - x(n-1))^2 + d(end)*(xq(i)-x(n-1))^3;
    i = i + 1;
end
    
yq(idx) = yq;

end

