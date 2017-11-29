function p = PolyFit( x,y,n )
%POLYFIT Summary of this function goes here
%   Detailed explanation goes here


x = x(:);
y = y(:);

% Construct Vandermonde matrix.
V(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
   V(:,j) = x.*V(:,j+1);
end

p = LS( V, y, 'QR');

end

