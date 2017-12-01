
n = 13;

f = @(x)(1./(1 + 25*x.^2));

xq = linspace(-1, 1, 1000);
y_true = f(xq);
plot(xq, y_true);
hold on

x = linspace(-1,1,n);
y = f(x);
yq = PolyInterp(x,y,xq);
plot(xq, yq);

yq = ChebyshevInterp(f, n, -1, 1, xq);
plot(xq, yq);

hold off

h = legend('$f(x)=\frac{1}{1+25x^2}$','equidistance interpolation','Chebyshev interpolation','location','south');
set(h,'Interpreter','latex');