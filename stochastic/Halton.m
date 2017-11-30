function u = Halton( p, n )
%HALTON Summary of this function goes here
%   Detailed explanation goes here

b = zeros(ceil(log(n) / log(p)), 1);
u = zeros(n, 1);
for jj = 1:n
    ii = 1;
    b(1) = b(1) + 1;
    while b(ii) > p - 1 + eps
        b(ii) = 0;
        ii = ii + 1;
        b(ii) = b(ii) + 1;
    end
    u(jj) = 0;
    for k = 1 : length(b)
        u(jj) = u(jj) + b(k) * p^(-k);
    end
end

end

