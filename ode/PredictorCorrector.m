function [x, y] = PredictorCorrector( f, a, b, y0, tol )
%PREDICTORCORRECTOR Summary of this function goes here
%   Detailed explanation goes here

hmax = (b - a) / 10;
h = hmax;
hmin = 1e-5;
x = zeros(11,1);
y = zeros(11,length(y0));
x(1) = a;
y(1, :) = y0;
flag = 1;
last = 0;

[xRK, yRK] = RK4(f, h, a, y0);
x(1:4) = xRK;
y(1:4, :) = yRK;
nflag = 1;
i = 5;
xnew = x(i-1) + h;

while(flag)
    if(i+3 > length(x))
        x = [x; zeros(length(x), 1)];
        y = [y; zeros(size(y))];
    end
    
    f3 = f(x(i-1), y(i-1,:));
    f2 = f(x(i-2), y(i-2,:));
    f1 = f(x(i-3), y(i-3,:));
    f0 = f(x(i-4), y(i-4,:));
    ypredict = y(i-1,:) + h*(55*f3 - 59*f2 + 37*f1 -9*f0) / 24;
    ycorrect = y(i-1,:) + h*(9*f(xnew, ypredict) + 19*f3 -5*f2 + f1)/24;
    sigma = 19*max(abs(ypredict - ycorrect)) / (270*h);
    if(sigma <= tol)
        y(i,:) = ycorrect;
        x(i) = xnew;
        if(last)
            flag = 0;
        else
            i = i + 1;
            nflag = 0;
            if(sigma <= 0.1*tol || x(i-1) + h > b)
                q = (tol / (2*sigma))^(1/4);
                if(q > 4)
                    h = 4*h;
                else
                    h = q*h;
                end
                if(h > hmax)
                    h = hmax;
                end
                if(x(i-1) + 4*h > b)
                    h = (b - x(i-1))/4;
                    last = 1;
                end
            
                [xRK, yRK] = RK4(f, h, x(i-1), y(i-1,:));
                x(i-1:i+2) = xRK;
                y(i-1:i+2,:) = yRK;
                nflag = 1;
                i = i + 3;
            end
        end
    else
        q = (tol / (2*sigma))^(1/4);
        if(q < 0.1 || isnan(q))
            h = 0.1*h;
        else
            h = q*h;
        end
        if(h < hmin)
            flag = 0;
            disp('hmin exceeded');
        else
            if(nflag)
                i = i - 3;
            end

            [xRK, yRK] = RK4(f, h, x(i-1), y(i-1,:));
            x(i-1:i+2) = xRK;
            y(i-1:i+2,:) = yRK;
            i = i + 3;
            nflag = 1;
        end
    end
    xnew = x(i-1) + h;
end

x = x(1:i);
y = y(1:i,:);

function [x, y] = RK4(f, h, x0, y0)
 
x = [x0; x0 + h; x0 + 2*h; x0 + 3*h];
y = zeros(4, length(y0));
y(1,:) = y0;

for i = 2:4
    k1 = h * f(x(i-1), y(i-1, :));
    k2 = h * f(x(i-1) + h/2, y(i-1,:) + k1/2);
    k3 = h * f(x(i-1) + h/2, y(i-1,:) + k2/2);
    k4 = h * f(x(i-1) + h, y(i-1,:) + k3);
    y(i,:) = y(i-1,:) + (k1 + 2*k2 + 2*k3 + k4)/6;
end

