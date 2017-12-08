function p = DogLeg( grad, delta, B )
%DOGLEG Summary of this function goes here
%   Detailed explanation goes here

delta2 = delta^2;
gg = grad'*grad;
gBg = grad'*B*grad;
pu = -gg/gBg * grad;
pu2 = pu'*pu;

if pu2 >= delta2
    p = delta/sqrt(pu2) * pu;
    return
end
pb = -SPDSolve(B, grad);
if Norm(pb) <= delta
    p = pb;
    return 
end
    
d = pb - pu;
a = d'*d;
b = pu'*d;
c = pu'*pu - delta2;

s = sqrt(b^2 - a*c);
r = (-b+s)/a;
if(r < 0 || r > 1)
    r = (-b-s)/a;
    if(r < 0 || r > 1)
        error(' ');
    end
end
p = pu + r*d;

    



end

