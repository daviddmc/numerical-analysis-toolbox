function [A, Q] = HessenbergReduction( A )
% QAQ^T = H;
%   Detailed explanation goes here

n = size(A, 1);
if(nargout > 1)
    Q = eye(n);
end
for ii = 1 : n-2
    v = House(A(ii+1:end, ii));
    A(ii+1:end, ii:end) = A(ii+1:end, ii:end) - 2*v*(v'*A(ii+1:end, ii:end));
    A(:,ii+1:end) = A(:,ii+1:end) - 2*(A(:,ii+1:end)*v)*v';
    if(nargout > 1)
        Q(ii+1:end, 1:end) = Q(ii+1:end, 1:end) - 2*v*(v'*Q(ii+1:end, 1:end));
    end
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


