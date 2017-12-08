function a = LineSearch( fun, grad, x, p, c1, c2, a1, aMax)
%LINESEARCH Summary of this function goes here
%   Detailed explanation goes here

% c1 = 1e-3
% c2 = 0.9
% a1 = 1;

a = a1;
i = 1;
f0 = fun(x);
g0 = grad(x)'*p;
fOld = f0;
gOld = g0;
aOld = 0;
while 1
    xNew = x+a*p;
    fNew = fun(xNew);
    if fNew > f0 + c1*a*g0 || (fNew >= fOld && i > 1)
        flag = 1;
        gNew = grad(xNew)'*p;
        break
    end
    gNew = grad(xNew)'*p;
    if abs(gNew) <= -c2*g0
        flag = 0;
        break;
    end
    if gNew >= 0
        flag = -1;
        break;
    end
    aOld = a;
    a = min(2*a, aMax);
    fOld = fNew;
    gOld = gNew;
    i = i+1;
end

if flag > 0
    aLow = aOld;
    fLow = fOld;
    gLow = gOld;
    aHigh = a;
    fHigh = fNew;
    gHigh = gNew;
elseif flag < 0
    aLow = a;
    fLow = fNew;
    gLow = gNew;
    aHigh = aOld;
    fHigh = fOld;
    gHigh = gOld;
end

if flag
    while 1
        d = aHigh - aLow;
        d1 = gLow + gHigh - 3*(fHigh - fLow)/(d);
        d2 = sign(d) * sqrt(d1^2 - gLow*gHigh);
        a = aHigh - d * (gHigh + d2 - d1) / (gHigh - gLow +2*d2);
        xNew = x+a*p;
        fNew = fun(xNew);
        if fNew > f0 + c1*a*g0 || fNew > fLow
            aHigh = a;
            fHigh = fNew;
            gHigh = grad(xNew)'*p;
        else
            gNew = grad(xNew)'*p;
            if abs(gNew) <= -c2*g0
                break;
            end
            if gNew*(aHigh - aLow) >= 0
                aHigh = aLow;
                fHigh = fLow;
                gHigh = gLow;
            end
            aLow = a;
            fLow = fNew;
            gLow = gNew;
        end
    end
end


end


