function [U,V,X,C,S] = GSVD(A,B,flag)

[m,p]  = size(A);
[n,pb] = size(B);
if pb ~= p
    error(message('MATLAB:gsvd:MatrixColMismatch'))
end

useQA = false;
useQB = false;
if nargin > 2
    if ~(isnumeric(flag) && isscalar(flag) && flag == 0)
        error(message('MATLAB:gsvd:InvalidFlag'));
    end
    % Economy-sized.
    useQA = m > p;
    useQB = n > p;
    if useQA
        [QA,A] = QR(A,0);
        m = p;
    end
    if useQB
        [QB,B] = QR(B,0);
        n = p;
    end
end

[Q,R] = QR([A;B],0);
[U,V,Z,C,S] = ThinCS(Q(1:m,:),Q(m+1:m+n,:));

if nargout < 2
    % Vector of generalized singular values.
    q = min(m+n,p);
    U = [zeros(q-m,1); diagk(C,max(0,q-m))]./[diagk(S,0); zeros(q-n,1)];
else
    % Full composition.
    X = R'*Z;
    if useQA
        U = QA*U;
    end
    if useQB
        V = QB*V;
    end
end

function D = diagk(X,k)
% DIAGK  K-th matrix diagonal.
% DIAGK(X,k) is the k-th diagonal of X, even if X is a vector.
if ~isvector(X)
    D = diag(X,k);
    D = D(:);  %Ensure column vector is returned for empty X.
else
    if ~isempty(X) && 0 <= k && 1+k <= size(X,2)
        D = X(1+k);
    elseif ~isempty(X) && k < 0 && 1-k <= size(X,1)
        D = X(1-k);
    else
        D = zeros(0,1,'like',X);
    end
end