function [A, Q, V] = BidiagonalReduction( A )
%

n = size(A, 1);
if(nargout > 1)
    Q = eye(n);
    if(nargout > 2)
        V = eye(n);
    end
end
for ii = 1 : n-2
    v = House(A(ii:end, ii));
    A(ii:end, ii:end) = A(ii:end, ii:end) - 2*v*(v'*A(ii:end, ii:end));
    if(nargout > 1)
        Q(ii:end, 1:end) = Q(ii:end, 1:end) - 2*v*(v'*Q(ii:end, 1:end));
    end
    
    v = House(A(ii, ii+1:end)');
    A(:,ii+1:end) = A(:,ii+1:end) - 2*(A(:,ii+1:end)*v)*v';
    if(nargout > 2)
        V(:,ii+1:end) = V(:,ii+1:end) - 2*(V(:,ii+1:end)*v)*v';
    end    
end

v = House(A(n-1:end, n-1));
A(n-1:end, n-1:end) = A(n-1:end, n-1:end) - 2*v*(v'*A(n-1:end, n-1:end));
if(nargout > 1)
    Q(n-1:end, 1:end) = Q(n-1:end, 1:end) - 2*v*(v'*Q(n-1:end, 1:end));
end


end

function x = House(x)
    if(x(1) == 0)
        x(1) = -Norm(x);
    else
        x(1) = x(1) + x(1)/abs(x(1)) * Norm(x);
    end
    x = x / Norm(x);
end


