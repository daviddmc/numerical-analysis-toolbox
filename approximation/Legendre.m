function P = Legendre( n )
%LEGENDRE Summary of this function goes here
%   Detailed explanation goes here

P0 = 1;
P1 = [0, 1];

if(n == 0)
    P = 1;
elseif(n==1)
    P = [0, 1];
else
    for ii = 2 : n
        P = [0  (2*ii - 1)*P1];
        P(1 : ii - 1) = P(1 : ii - 1) - (ii-1) * P0;
        P = P / ii;
        P0 = P1;
        P1 = P;
    end
       
end

