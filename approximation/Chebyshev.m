function T = Chebyshev( n )

T1 = 1;
T2 = [0, 1];

if(n == 0)
    T = T1;
elseif(n == 1)
    T = T2;
else 
    for ii = 2 : n
        T = [0  2*T2];
        T(1 : ii - 1) = T(1 : ii - 1) - T1; 
        T1 = T2;
        T2 = T;
    end
end

