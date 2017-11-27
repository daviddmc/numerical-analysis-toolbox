function [A, J] = Jacobi(A, epsilon)

n = size(A, 1);
delta = Norm(A, 'fro') / n;

if nargout>1
    J = eye(n);
end

while(delta > epsilon)
    for j = 1 : n
        for k = 2 : n
            if(abs(A(j,k)) > delta)
                tau = (A(j,j) - A(k,k)) / (2 * A(j,k));
                t = sign(tau) / (abs(tau) + sqrt(1 + tau^2));
                c = 1 / sqrt(1 + t^2);
                s = c*t;
                q = [c s; -s c];
                A([j, k], :) = q * A([j, k],:);
                A(:, [j, k]) = A(:, [j, k]) * q';
                if nargout > 1
                    J(:, [j, k]) = J(:, [j, k]) * q';
                end
            end
        end
    end
    delta = delta / n;
end




