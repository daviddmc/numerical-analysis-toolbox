function r = Brent( f, a, b )
%BRENT Summary of this function goes here
%   Detailed explanation goes here


delta = 1e-6;
fa = f(a);
fb = f(b);
r = [];
if(fa*fb > 0)
    return
end

if(abs(fa) < abs(fb))
    temp = a;
    a = b;
    b = temp;
    temp = fa;
    fa = fb;
    fb = temp;
end

c = a;
fc = f(c);
mflag = 1;

while(1)
    if(fa ~= fc && fb ~= fc)
        s = a * fb * fc / (fa - fb) / (fa - fc);
        s = s + b * fa * fc / (fb - fa) / (fb - fc);
        s = s + c * fa * fb / (fc - fa) / (fc - fb);
    else
        s = b - fb * (b - a) / (fb - fa);
    end
    
    if((s - (3*a + b) / 4) * (s - b) > 0 || ...
            (mflag && abs(s-b) >= abs(b-c)/2) ||...
            (~mflag && abs(s-b) >= abs(c-d)/2) ||...
            (mflag && abs(b-c) < delta) ||...
            (~mflag && abs(c-d) < delta))
        s = (a+b)/2;
        mflag = 1;
    else
        mflag = 0;
    end
    fs = f(s);
    d = c;
    c = b;
    if(fa*fs < 0)
        b = s;
        fb = f(b);
    else
        a = s;
        fa = f(a);
    end
    
    if(abs(fa) < abs(fb))
        temp = a;
        a = b;
        b = temp;
        temp = fa;
        fa = fb;
        fb = temp;
    end
    
    if(fb == 0)
        r = b;
        break;
    elseif(fs == 0 || abs(b - a) < delta)
        r = s;
        break
    end
    
end



end

