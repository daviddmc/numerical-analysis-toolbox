function I = Romberg( f, a, b, kmax, epsilon)
%TRAPEZOID Summary of this function goes here
%   Detailed explanation goes here

T = zeros(kmax+1, 1);
T(1) = Trapezoid(f, a, b, 1);
I = T(1);
for k = 1 : kmax
    T1 = Trapezoid(f, a, b, 2^k);
    for jj = 1 : k+1
        T2 = (4^jj * T1 - T(jj)) / (4^jj - 1);
        T(jj) = T1;
        T1 = T2;
    end
    if(abs(I - T(k+1)) < epsilon)
        I = T(k+1);
        break
    end
    I = T(k+1);
end

end

