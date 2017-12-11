function [Q, R, E] = QR(A, varargin)
% QR     Orthogonal-triangular decomposition.
%   [Q,R] = QR(A), where A is m-by-n, produces an m-by-n upper triangular
%   matrix R and an m-by-m unitary matrix Q so that A = Q*R.
%
%   [Q,R] = QR(A,0) produces the "economy size" decomposition.
%   If m>n, only the first n columns of Q and the first n rows of R are
%   computed. If m<=n, this is the same as [Q,R] = qr(A).
%
%   [Q,R,E] = QR(A) produces unitary Q, upper triangular R and a
%   permutation matrix E so that A*E = Q*R. The column permutation E is
%   chosen so that ABS(DIAG(R)) is decreasing.
%   
%   [Q,R,e] = QR(A,'vector') returns the permutation information as a
%   vector instead of a matrix.  That is, e is a row vector such that 
%   A(:,e) = Q*R. Similarly, [Q,R,E] = qr(A,'matrix') returns a permutation 
%   matrix E. This is the default behavior.
% 
%   [Q,R,e] = QR(A,0) produces an "economy size" decomposition in which e
%   is a permutation vector, so that A(:,e) = Q*R.
%
%   [Q,R,...] = QR(..., METHOD) specifies the algorithm for QR
%   decomposition. The available method are:
%
%     'Householder' - (default) Householder reflection
%     'Givens'      - Givens rotation
%     'GramSchmidt' - Gram-Schmidt orthogonalization, this method will
%                     always produce "economy size" decomposition. 
%
%   See also LU, Cholesky.

%   Copyright 2017 Junshen Xu
    
flagPivot = (nargout == 3);
[method, flagVec, flagEcon] = ParseInputString(varargin);
if flagPivot && flagEcon
    flagVec = 1;
end

if(strncmpi(method, 'GramSchmidt',2))
    [r, c] = size(A);
    
    if flagPivot
        p = 1:c;
    end
    Q = A;
    if(r >= c)
        R = zeros(c);
    else
        R = zeros(size(A));
    end
    for k = 1 : min(r,c)
        if flagPivot
            C = sum(Q(:,k:end).^2);
            [maxNormSquare, pivotK] = max(C);
            pivotK = pivotK + k - 1;
            p([pivotK, k]) = p([k, pivotK]);
            Q(:, [pivotK, k]) = Q(:, [k, pivotK]);
            R(k, k) = sqrt(maxNormSquare);
            if(R(k, k) == 0)
                break;
            end
        else
            R(k, k) = Norm(Q(:, k));
            if(R(k, k) == 0)
                continue;
            end
        end
        Q(:, k) = Q(:, k) / R(k, k);
        R(k, k+1:end) = Q(:,k)'*Q(:,k+1:end);
        if flagPivot
            R(1:k-1, [pivotK, k]) = R(1:k-1, [k, pivotK]);
        end
        Q(:,k+1:end) = Q(:,k+1:end) - Q(:,k)*R(k,k+1:end);
    end
    if r < c
        Q = Q(:, 1:r);
    end
elseif(strncmpi(method, 'Givens',2))
    [m, n] = size(A);
    if flagPivot
        C = sum(abs(A).^2);
        p = 1:n;
    end
    R = A;
    for ii = 1 : n
        if flagPivot
            [tau, pivotIdx] = max(C(ii:end));
            pivotIdx = pivotIdx + ii - 1;
            if(tau <= 0)
                break;
            end
            p([pivotIdx, ii]) = p([ii, pivotIdx]);
            C([pivotIdx, ii]) = C([ii, pivotIdx]);
            R(:, [pivotIdx, ii]) = R(:, [ii, pivotIdx]);
        end
        for k = ii+1 : m
            if(R(k, ii) == 0)
                continue;
            end
            t = sqrt(abs(R(ii, ii))^2 + abs(R(k, ii))^2);
            c = abs(R(ii,ii)) / t;
            if c == 0
                s = R(k, ii)' / t;
            else
                s = R(k, ii)' * sign(R(ii,ii)) / t;
            end
            R([ii, k], ii:end) = [c, s; -s', c] * R([ii, k], ii:end);
            R(k, ii) = s;
        end
        if flagPivot
            C(ii+1:end) = C(ii+1:end) - abs(R(ii,ii+1:end)).^2;
        end
    end
    if flagEcon && m > n
        Q = eye(m, n);
    else
        Q = eye(m);
    end
    for ii = ii : -1 : 1
        for k = m : -1 : ii+1 
            rho = R(k, ii);
            if rho == 0
                continue;
            end
            s = rho;
            c = sqrt(1 - abs(rho)^2);
            Q([ii, k], ii:end) = [c, -s; s', c] * Q([ii, k], ii:end);  
        end
    end
    R = triu(R);
    if flagEcon && m > n
        R = R(1:n,:);
    end
elseif(strncmpi(method,'Householder',1))
    [r, c] = size(A);
    if flagPivot
        C = sum(abs(A).^2);
        p = 1:c;
    end
    beta = zeros(c, 1);
    for ii = 1 : min(r,c)
        if flagPivot
            [tau, pivotIdx] = max(C(ii:end));
            pivotIdx = pivotIdx + ii - 1;
            if(tau <= 0)
                break;
            end
            p([pivotIdx, ii]) = p([ii, pivotIdx]);
            C([pivotIdx, ii]) = C([ii, pivotIdx]);
            A(:, [pivotIdx, ii]) = A(:, [ii, pivotIdx]);
        end
        [v, beta(ii)] = HouseVec(A(ii:end, ii));
        A(ii:end, ii:end) = A(ii:end, ii:end) - beta(ii)*v*(v'*A(ii:end, ii:end));
        A(ii+1:end, ii) = v(2:end);  
        if flagPivot
            C(ii+1:end) = C(ii+1:end) - abs(A(ii,ii+1:end)).^2;
        end
    end
    if flagEcon && r > c
        Q = eye(r, c);
    else
        Q = eye(r);
    end
    for ii = ii : -1 : 1
        v = [1; A(ii+1:end, ii)];
        Q(ii:end, ii:end) = Q(ii:end, ii:end) - beta(ii)*v*(v'*Q(ii:end, ii:end));
    end
    R = triu(A);
    if flagEcon && r > c
        R = R(1:c,:);
    end
else
    error('method error');
end

if flagPivot
    if flagVec
        E = p;
    else
        E = eye(size(A, 2));
        E = E(:, p);
    end
end

function [method, flagVec, flagEcon] = ParseInputString(stringCell)

method = 'h';
flagVec = 0;
flagEcon = 0;
 
if length(stringCell) == 1
    a = stringCell{1};
    if a == 0
        flagEcon = 1;
        return
    end
    if(~ischar(a)&&~isstring(a))
        error('extra input of QR should be 0 or string');
    end
    if(strncmpi(a, 'vector', 1))
        flagVec = 1;
        return
    end
    if(~strncmpi(a, 'matrix', 1))
        method = a;
        return
    end
elseif length(stringCell) == 2
    if stringCell{1} == 0
        flagEcon = 1;
    elseif strncmpi(stringCell{1}, 'vector', 1)
        flagVec = 1;
    end
    method = stringCell{2};
elseif length(stringCell) > 2
    error('too many input');
end

function [x, beta] = HouseVec(x)

sigma = Norm(x);
if(sigma == 0)
    beta = 0;
    return;
else
    beta = 1 / (sigma*(sigma + abs(x(1))));
end
    
if x(1) == 0
    mu = sigma;
else
    mu = -x(1)/abs(x(1)) * sigma;
end

x(1) = x(1) - mu;
beta = beta * abs(x(1))^2;
x = x / x(1);

    

