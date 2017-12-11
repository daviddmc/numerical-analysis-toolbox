function X = WeightedLS( A, B, D )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


D = sqrt(D);
A = D.*A;
B = D.*B;

X = LS(A, B);

end

