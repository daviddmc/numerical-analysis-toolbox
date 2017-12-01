function [xa, ya] = RungeKuttaFehlberg( f, a, b, y0, tol )
%RUNGEKUTTAFEHLBERG Summary of this function goes here
%   Detailed explanation goes here

xa = a;
ya = y0;

x = a;
y = y0;
hmax = (b - a) / 10;
hmin = 1e-6;
h = hmax;
flag = 1;

while(flag)
    k1 = h * f(x, y);
    k2 = h * f(x + h/4, y + k1 / 4);
    k3 = h * f(x + 3*h/8, y + 3*k1/32 + 9*k2/32);
    k4 = h * f(x + 12*h/13, y + 1932*k1/2197 - 7200*k2/2197 + 7296*k3/2197);
    k5 = h * f(x + h, y + 439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104);
    k6 = h * f(x + h/2, y - 8*k1/27 + 2*k2 - 3544*k3/2565 + 1859*k4/4104 - 11*k5/40);
    R = abs(k1/360 - 128*k3/4275 - 2197*k4/75240 + k5/50 + 2*k6/55) / h;
    
    if(R <= tol)
        x = x+h;
        xa = [xa x];
        y = y + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5;
        ya = [ya, y];
    end
    
    delta = 0.84 * (tol / R)^(1/4);
    if(delta <= 0.1)
        h = h / 10;
    elseif(delta >= 4)
        h = 4 * h;
    else
        h = delta * h;
    end
    
    if(h > hmax)
        h = hmax;
    end
    
    if(x >= b)
        flag = 0;
    elseif(x + h > b)
        h = b - x;
    elseif(h < hmin)
        flag = 0;
        disp('minimum h exceeded');
    end
    


end

