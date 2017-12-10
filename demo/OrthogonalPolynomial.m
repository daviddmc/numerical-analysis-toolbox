% This script plots some common orthogonal polynomials.
%
% Copyright 2017 Junshen Xu

%% setup
N = 4;

%% Chebyshev
figure
hold on
x = linspace(-1, 1, 100);
label = {};
for i = 0 : N
    p = Chebyshev(i);
    y = Horner(p, x);
    plot(x, y);
    label = [label sprintf('$n=%d$', i)];
end
hold off
title('Chebyshev polynomials of first type', 'interpreter', 'latex');
xlabel('$x$','interpreter', 'latex');
ylabel('$y$','interpreter', 'latex');
h = legend(label,'location','south');
set(h,'Interpreter','latex');
grid on

%% Legendre
figure
hold on
x = linspace(-1, 1, 100);
label = {};
for i = 0 : N
    p = Legendre(i);
    y = Horner(p, x);
    plot(x, y);
    label = [label sprintf('$n=%d$', i)];
end
hold off
title('Legendre polynomials', 'interpreter', 'latex');
xlabel('$x$','interpreter', 'latex');
ylabel('$y$','interpreter', 'latex');
h = legend(label,'location','south');
set(h,'Interpreter','latex');
grid on

%% Laguerre
figure
hold on
x = linspace(0, 5, 100);
label = {};
for i = 0 : N
    p = Laguerre(i);
    y = Horner(p, x);
    plot(x, y);
    label = [label sprintf('$n=%d$', i)];
end
hold off
title('Laguerre polynomials', 'interpreter', 'latex');
xlabel('$x$','interpreter', 'latex');
ylabel('$y$','interpreter', 'latex');
h = legend(label,'location','south');
set(h,'Interpreter','latex');
grid on

%% Hermite
figure
hold on
x = linspace(-2, 2, 100);
label = {};
for i = 0 : N
    p = Hermite(i);
    y = Horner(p, x);
    plot(x, y);
    label = [label sprintf('$n=%d$', i)];
end
hold off
title('Hermite polynomials', 'interpreter', 'latex');
xlabel('$x$','interpreter', 'latex');
ylabel('$y$','interpreter', 'latex');
h = legend(label,'location','north');
set(h,'Interpreter','latex');
grid on

