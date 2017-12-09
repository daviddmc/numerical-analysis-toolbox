% This script shows the Runge phenomenon Runge's phenomenon which is a 
% problem of oscillation at the edges of an interval that occurs when using
% polynomial interpolation with polynomials of high degree over a set of 
% equispaced interpolation points.

% Copyright 2017 Junshen Xu


ns = [11, 13, 15]; % orders of polynomials
f = @(x)(1./(1 + 25*x.^2)); % the Runge function.
% plot the Runge funciton
xq = linspace(-1, 1, 1000);
y_true = f(xq);
plot(xq, y_true);
label = {'$f(x)=\frac{1}{1+25x^2}$'};
% plot interpolation results
hold on
for n = ns
    x = linspace(-1,1,n);
    y = f(x);
    yq = PolyInterp(x,y,xq);
    plot(xq, yq);
    label = [label sprintf('$n=%d$', n)];
end
hold off
title('Runge effect')
xlabel('x');
ylabel('y');
h = legend(label,'location','best');
set(h,'Interpreter','latex');