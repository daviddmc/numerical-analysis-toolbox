function r = PolyRealRoot( a )
%POLYREALROOT Summary of this function goes here
%   Detailed explanation goes here

n = length(a) - 1;
a = a(end:-1:1);
p = a;
r = [];
z = 0;

for k = n : -1 : 1
    pp = p(1 : k) .* (k : -1: 1);
    
    [z, flag] = Newton( @(x)(PolyValue(p , x)), ...
        @(x)(PolyValue(pp, x)), ...
        z, 20, 1e-6);
    
    p = filter(1,[1 -z],p);
    p(end) = [];
    
    if(flag)
        r = [r z];
    else
        break;
    end
end

p = a;
pp = p(1 : n) .* (n : -1: 1);
f = @(x)(PolyValue(p , x));
df = @(x)(PolyValue(pp, x));
for ii = 1 : length(r)
    r(ii) = Newton(f,df,r(ii),5, 1e-10);
end

end

function y = PolyValue(p, x)
    y = filter(1,[1 -x],p);
    y = y(length(p));
end


