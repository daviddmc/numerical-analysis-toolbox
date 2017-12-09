% This script shows the round-off error in numerical differentiation

% step size
h = logspace(0, -15,100);
f0 = cos(pi/6);
fh = cos(pi/6 + h);
% differentiation
dfdx = (fh - f0) ./ h;
% error
err = abs(dfdx + 0.5);
% plot results
loglog(h, err);
title('round-off error in numerical differentiation')
set(gca,'xdir','reverse')
xlabel('step size')
ylabel('error')
grid on
