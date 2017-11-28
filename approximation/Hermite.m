function H = Hermite( n)
%HERMITE Summary of this function goes here
%   Detailed explanation goes here
H0 = 1;
H1 = [0, 2];

if(n == 0)
    H = H0;
elseif(n == 1)
    H = H1;
else 
    for ii = 2 : n
        H = [0 2*H1];
        H(1 : ii - 1) = H(1 : ii - 1) - (ii-1)*2 * H0; 
        H0 = H1;
        H1 = H;
    end
end