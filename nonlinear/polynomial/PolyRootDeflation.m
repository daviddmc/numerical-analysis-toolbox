function r = PolyRootDeflation( a )
% PolyRootDeflation    Find roots of polynomial using deflation.
%   PolyRootDeflation(A) computes the roots of the polynomial using Muller 
%   method and deflation. The coefficients vector A has N+1 components, 
%   the polynomial is A(1)*X^N + ... + A(N)*X + A(N+1).
%
%   See also PolyRoot, PolyRealRootDeflation, Muller.

%   Copyright 2017 Junshen Xu

a = a(:).';
firstNonZero = find(a, 1);
lastNonZero = find(a, 1, 'last');
if isempty(firstNonZero)
    r = [];
    return
end
a = a(firstNonZero:lastNonZero);
r = zeros(1, lastNonZero - firstNonZero - 1);

d = a(2:end)./a(1);
while any(isinf(d))
    a = a(2:end);
    d = a(2:end)./a(1);
end

n = length(a) - 1;
p = a;
rs = [];
z = 0;

epsilon1 = 1e-6;
epsilon2 = 1e-10;

flagReal = isreal(a);

k = n;
while(k > 0)
    [z, flag]  = Muller(@(x)(PolyValue(p , x)), ...
        z+100*epsilon1, z-100*epsilon1, z, epsilon1, 100);
    if(flag == 0)
        if flagReal
            if(abs(imag(z)) < epsilon1)
                rs = [rs real(z)];
                p = filter(1,[1 -z],p);
                p(end) = [];
                k = k - 1;
            else
                rs = [rs z z'];
                p = filter(1,[1 -z],p);
                p(end) = [];
                p = filter(1,[1 -z'],p);
                p(end) = [];
                p = real(p);
                k = k - 2;
            end
        else
            rs = [rs z];
            p = filter(1,[1 -z],p);
            p(end) = [];
            k = k - 1;
        end
    else
        break;
    end        
end

p = a;
f = @(x)(PolyValue(p , x));
k = length(rs);
while(k > 0)
    if(isreal(rs(k)) || ~flagReal)
        rs(k) = Muller(f, rs(k) + 100*epsilon2, rs(k) - 100*epsilon2, rs(k),...
            epsilon2, 5);
        k = k - 1;
    else
        rs(k) = Muller(f, rs(k) + 100*epsilon2, rs(k) - 100*epsilon2, rs(k),...
            epsilon2, 5);
        rs(k - 1) = rs(k)';
        k = k - 2;
    end
end

r = [r rs];

end

function y = PolyValue(p, x)
    y = filter(1,[1 -x],p);
    y = y(length(p));
end


