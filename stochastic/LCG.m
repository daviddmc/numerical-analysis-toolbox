function x = LCG( n, a, b, m, seed, interval )
%LCG Summary of this function goes here
%   Detailed explanation goes here

m = 2147483647; %2^31 -1 
a = 16807;
b = 0;
seed = mod(round(100 * cputime), m);
interval = [0, 1];

x = zeros(n,1);
x(1) = mod(a*seed+b,m);

for ii = 2 : n
    x(ii) = mod(a * x(ii-1) + b, m);
end

x = x / m;
x = interval(1) + (interval(2) - interval(1)) * x;

end

