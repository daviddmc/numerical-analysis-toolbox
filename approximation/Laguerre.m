function L = Laguerre( n )
%LAGUERRE Summary of this function goes here
%   Detailed explanation goes here
L0 = 1;
L1 = [1, -1];

if(n == 0)
    L = L0;
elseif(n == 1)
    L = L1;
else 
    for ii = 2 : n
        L = [(2*ii-1) * L1 0] - [0 L1];
        L(1 : ii - 1) = L(1 : ii - 1) - (ii-1)^2 * L0; 
        L0 = L1;
        L1 = L;
    end
end


