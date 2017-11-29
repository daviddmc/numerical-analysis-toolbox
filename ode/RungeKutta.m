function [x, y] = RungeKutta( f, a, b, y0 )
%RUNGEKUTTA Summary of this function goes here
%   Detailed explanation goes here

N = 10;
hmax = (b - a) / N;
epsilon = 1e-7;
x = zeros(N+1, 1);
x(1) = a;
y = zeros(N+1, 1);
y(1) = y0;
n = 0;
wt = max(abs(y0),1e-3);


h = hmax;
xdel = (a + min(sqrt(eps)*max(abs(a),abs(a+h)),h)) - a;
f1 = f(a+xdel, y0); 
f0 = f(a, y0);
dfdt = (f1 - f0) ./ xdel;
yp = f0;
f1 = f(a, y0+1e-6);
dfdy = (f1 - f0) / (1e-6);
DfDt = dfdt + dfdy*yp;
rh = 1.43 * sqrt(0.5 * (norm(DfDt) / wt) / 1e-3);

if h * rh > 1
    h = 1 / rh;
end

h = min(h, hmax);

while(1)
    for n = (n+1) : N
        if(x(n) + h > b)
            h = b - x(n);
        end
        
        if(x(n) >= b)
            break;
        end

        flag = 0;

        while(1)
            if(flag == -1)% half
                yh = yhalf;
                yhalf = Update(f, x(n), y(n), h/2);
            elseif(flag == 1)% doubel
                yhalf = yh;
                yh = Update(f, x(n), y(n), h);
            elseif(flag == 0)% init
                yh = Update(f, x(n), y(n), h);
                yhalf = Update(f, x(n), y(n), h/2);
            end
            yhalf2 = Update(f, x(n) + h/2, yhalf, h/2);

            if(abs(yh - yhalf2) < 10 * epsilon)
                if(flag == -1)
                    x(n+1) = x(n) + h;
                    y(n+1) = yh;
                    break;
                else
                    if(h * 2 > hmax || flag == 1)
                        x(n+1) = x(n) + h;
                        y(n+1) = yh;
                        break;
                    end
                    
                    h = 2*h;
                    
                    if(x(n) + h > b)
                        x(n+1) = b;
                        y(n+1) = Update(f, x(n), y(n), b - x(n));
                        break;
                    end
                    flag = 1;
                end
            else
                h = h/2;
                if(flag == 1)
                    x(n+1) = x(n) + h;
                    y(n+1) = yhalf;
                    break;
                else
                    flag = -1;
                end
            end
        end
    end
    
    if(n < N || x(n) >= b)
        x = x(1:n);
        y = y(1:n);
        break;
    end

    x = [x ; zeros(N, 1)];
    y = [y ; zeros(N, 1)];
    N = 2*N;
end
    
    
end


function y = Update(f, xn, yn, h)
    k1 = f(xn, yn);
    k2 = f(xn+h/2, yn + h*k1/2);
    k3 = f(xn+h/2, yn + h*k2/2);
    k4 = f(xn+h, yn+h*k3);
    y = yn + h*(k1 + 2*k2 + 2*k3 + k4)/6;
end

