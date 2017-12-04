function p = PolyFit( x,y,n )
%PolyFit    Least squares fit polynomial to data.
%   P = PolyFit(X, Y, N) finds the coefficients of a polynomial P of degree
%   N that fits the data Y best in a least-squares sense. P is a row vector
%   of length N+1 containing the polynomial coefficients in desending
%   powers, P(1)*X^N + P(2)*X^(N-1) + ... + P(N)*X + P(N+1).

x = x(:);
y = y(:);

% Construct Vandermonde matrix.
V(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
   V(:,j) = x.*V(:,j+1);
end

% solve ls problem
p = LS( V, y, 'QR');

end

