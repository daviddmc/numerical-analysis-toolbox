function yq = Neville( x, y, Xq)
%NEVILLE Summary of this function goes here
%   Detailed explanation goes here

%epsilon = 1e-6;
n = length(x);
Q = zeros(n);
Q(:, 1) = y;
yq = zeros(size(Xq));
for k = 1 : length(yq)
    xq = Xq(k);
    for i = 1 : n - 1
        for j = 1 : i
            Q(i+1, j+1) = ((xq - x(i-j+1)) * Q(i+1, j) - (xq - x(i+1)) * Q(i, j)) / (x(i+1) - x(i-j+1));
        end
    end
    yq(k) = Q(n ,n);
end


end

