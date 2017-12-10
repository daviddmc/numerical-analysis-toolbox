function r = PolyRootDeflate( a )
%POLYREALROOT Summary of this function goes here
%   Detailed explanation goes here

n = length(a) - 1;
%a = a(end:-1:1);
p = a;
r = [];
z = 0;

epsilon1 = 1e-6;
epsilon2 = 1e-10;

k = n;
while(k > 0)
    [z, flag] = Muller(@(x)(PolyValue(p , x)), ...
        z, 100, epsilon1);
    if(flag)
        if(abs(imag(z)) < epsilon1)
            r = [r real(z)];
            p = filter(1,[1 -z],p);
            p(end) = [];
            k = k - 1;
        else
            r = [r z z'];
            p = filter(1,[1 -z],p);
            p(end) = [];
            p = filter(1,[1 -z'],p);
            p(end) = [];
            p = real(p);
            k = k - 2;
        end
    else
        break;
    end
    
    if(isreal(z))
        
    else
        
    end
        
end

p = a;
f = @(x)(PolyValue(p , x));
k = length(r);
while(k > 0)
    if(isreal(r(k)))
        r(k) = Muller(f, r(k), 5, epsilon2);
        k = k - 1;
    else
        r(k) = Muller(f, r(k), 5, epsilon2);
        r(k - 1) = r(k)';
        k = k - 2;
    end
end

end

function y = PolyValue(p, x)
    y = filter(1,[1 -x],p);
    y = y(length(p));
end


