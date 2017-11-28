
ns = [11, 13, 15];

f = @(x)(1./(1 + 25*x.^2));

xq = linspace(-1, 1, 1000);
y_true = f(xq);
plot(xq, y_true);
hold on

label = {'$f(x)=\frac{1}{1+25x^2}$'};
for n = ns
    x = linspace(-1,1,n);
    y = f(x);
    yq = LagrangeInterp(x,y,xq);
    plot(xq, yq);
    label = [label sprintf('$n=%d$', n)];
end

hold off

h = legend(label,'location','best');
set(h,'Interpreter','latex');