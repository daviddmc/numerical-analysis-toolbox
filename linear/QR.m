function [Q, R] = QR(A, method)
% QR decomposition
% input
% A
% method: 'GramSchmidt', 'Givens' or 'Householder', default: Householder
% output
% Q
% R

if(~exist('method','var'))
    method = 'Householder';
end

if(strcmp(method, 'GramSchmidt'))
    r = size(A, 1);
    c = size(A, 2);
    if(c < r)
        %error(' numRow must be not more than numCol ')
    end
    
    if(r >= c)
        Q = zeros(size(A));
        R = zeros(c);
        for ii = 1 : c
            u = A(:, ii);
            for jj = 1 : ii - 1
                R(jj, ii) = Q(:, jj)' * u; 
                u = u - R(jj, ii) * Q(:, jj);
            end
            R(ii, ii) = Norm(u);
            Q(:, ii) = u / R(ii, ii);
        end
    else
        Q = zeros(r);
        R = zeros(size(A));
        for ii = 1 : r
            u = A(:, ii);
            for jj = 1 : ii - 1
                R(jj, ii) = Q(:, jj)' * u; 
                u = u - R(jj, ii) * Q(:, jj);
            end
            R(ii, ii) = Norm(u);
            Q(:, ii) = u / R(ii, ii);
        end
        R(:, r+1:c) = Q' * A(:, r+1:c);
    end

    %R = Q' * A;
    
elseif(strcmp(method, 'Givens'))
    r = size(A, 1);
    Q = eye(r);
    R = A;
    for ii = 1 : r-1
        for k = ii+1 : r
           t = sqrt(R(ii, ii)^2 + R(k, ii)^2);
           if(t == 0)
               continue
           end
           c = R(ii,ii) / t;
           s = R(k, ii) / t;
           R([ii, k], :) = [c, s; -s, c] * R([ii, k], :);
           Q([ii, k], :) = [c, s; -s, c] * Q([ii, k], :);
           %Q(:, [ii, k]) = Q(:, [ii, k]) * [c, -s; s, c]; 
        end
    end
    Q = Q';
    
elseif(strcmp(method,'Householder'))
    r = size(A, 1);
    Q = eye(r);
    R = A;
    for ii = 1 : r-1
        x = R(ii:end, ii);
        e = zeros(size(x));
        if(x(1) == 0)
            e(1) = Norm(x);
        else
            e(1) = -x(1)/abs(x(1)) * Norm(x);
        end
        u = x - e;
        v = u / Norm(u);
        R(ii:end, ii:end) = R(ii:end, ii:end) - 2*v*(v'*R(ii:end, ii:end));
        Q(:,ii:end) = Q(:,ii:end) - 2*(Q(:,ii:end)*v)*v';     
    end
else
    error('method error');
end

