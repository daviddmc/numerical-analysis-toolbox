function r = PolyRealRootDeflation( a )
% PolyRealRootDeflation    Find real roots of real-coefficient polynomial.
%   roots(A) computes the real roots of the polynomial using Newton method
%   and deflation. The coefficients vector A should by real. If A has N+1 
%   components, the polynomial is A(1)*X^N + ... + A(N)*X + A(N+1).
%
%   See also PolyRoot, PolyRootDeflation, Newton.

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

for k = n : -1 : 1
    pp = PolyDiff(p);
    [z, flag] = Newton( @(x)(PolyValue(p , x)), ...
        @(x)(PolyValue(pp, x)), ...
        z, 1e-8, 30);
    p = filter(1,[1 -z],p);
    p(end) = [];
    if(flag == 0)
        rs = [rs z];
    else
        warning('The polynomial may have complex roots.');
        break;
    end
end

p = a;
pp = PolyDiff(p);
f = @(x)(PolyValue(p , x));
df = @(x)(PolyValue(pp, x));
for ii = 1 : length(rs)
    rs(ii) = Newton(f,df,rs(ii),1e-10,5);
end

r = [r rs];

end

function y = PolyValue(p, x)
    y = filter(1,[1 -x],p);
    y = y(length(p));
end


