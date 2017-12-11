function [U1,V1,V,C,S] = ThinCS(Q1,Q2)
% CSD  Cosine-Sine Decomposition
% [U,V,Z,C,S] = csd(Q1,Q2)
%
% Given Q1 and Q2 such that Q1'*Q1 + Q2'*Q2 = I, the
% C-S Decomposition is a joint factorization of the form
%    Q1 = U*C*Z' and Q2=V*S*Z'
% where U, V, and Z are orthogonal matrices and C and S
% are diagonal matrices (not necessarily square) satisfying
%    C'*C + S'*S = I

[m,p] = size(Q1);
n = size(Q2,1);

if m < n
    [V1,U1,V,S,C] = ThinCS(Q2,Q1);
    j = p:-1:1; C = C(:,j); S = S(:,j); V = V(:,j);
    m = min(m,p); i = m:-1:1; C(1:m,:) = C(i,:); U1(:,1:m) = U1(:,i);
    n = min(n,p); i = n:-1:1; S(1:n,:) = S(i,:); V1(:,1:n) = V1(:,i);
    return
end
% Henceforth, n <= m.

[C,U1,V] = SVD(Q1);

q = min(m,p);
i = 1:q;
j = q:-1:1;
C(i,i) = C(j,j);
U1(:,i) = U1(:,j);
V(:,i) = V(:,j);
S = Q2*V;

if q == 1
    k = 0;
elseif m < p
    k = n;
else
    k = max([0; find(diag(C) <= 1/sqrt(2))]);
end
[V1,~] = QR(S(:,1:k));
S = V1'*S;
r = min(k,m);
S(:,1:r) = diagf(S(:,1:r));
if m == 1 && p > 1
    S(1,1) = 0;
end

if k < min(n,p)
    r = min(n,p);
    i = k+1:n;
    j = k+1:r;
    [ST,UT,VT] = SVD(S(i,j));
    if k > 0
        S(1:k,j) = 0;
    end
    S(i,j) = ST;
    C(:,j) = C(:,j)*VT;
    V1(:,i) = V1(:,i)*UT;
    V(:,j) = V(:,j)*VT;
    i = k+1:q;
    [Q,R] = QR(C(i,j));
    C(i,j) = diagf(R);
    U1(:,i) = U1(:,i)*Q;
end

if m < p
    % Diagonalize final block of S and permute blocks.
    q = min([nnz(abs(diagk(C,0))>10*m*eps(class(C))), ...
        nnz(abs(diagk(S,0))>10*n*eps(class(C))), ...
        nnz(max(abs(S(:,m+1:p)),[],2)<sqrt(eps(class(C))))]);
    
    % maxq: maximum size of q such that the expression used later on,
    %        i = [q+1:q+p-m, 1:q, q+p-m+1:n],
    % is still a valid permutation.
    maxq = m+n-p;
    q = q + nnz(max(abs(S(:,q+1:maxq)),[],1)>sqrt(eps(class(C))));
    
    i = q+1:n;
    j = m+1:p;
    % At this point, S(i,j) should have orthogonal columns and the
    % elements of S(:,q+1:p) outside of S(i,j) should be negligible.
    [Q,R] = QR(S(i,j));
    S(:,q+1:p) = 0;
    S(i,j) = diagf(R);
    V1(:,i) = V1(:,i)*Q;
    if n > 1
        i = [q+1:q+p-m, 1:q, q+p-m+1:n];
    else
        i = 1;
    end
    j = [m+1:p 1:m];
    C = C(:,j);
    S = S(i,j);
    V = V(:,j);
    V1 = V1(:,i);
end

if n < p
    % Final block of S is negligible.
    S(:,n+1:p) = 0;
end

% Make sure C and S are real and positive.
[U1,C] = diagp(U1,C,max(0,p-m));
C = real(C);
[V1,S] = diagp(V1,S,0);
S = real(S);

% ------------------------

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

% ------------------------

function X = diagf(X)
% DIAGF  Diagonal force.
% X = DIAGF(X) zeros all the elements off the main diagonal of X.
X = triu(tril(X));

% ------------------------

function [Y,X] = diagp(Y,X,k)
% DIAGP  Diagonal positive.
% [Y,X] = diagp(Y,X,k) scales the columns of Y and the rows of X by
% unimodular factors to make the k-th diagonal of X real and positive.
D = diagk(X,k);
j = find(real(D) < 0 | imag(D) ~= 0);
D = diag(conj(D(j))./abs(D(j)));
Y(:,j) = Y(:,j)*D';
X(j,:) = D*X(j,:);
X = X+0; % use "+0" to set possible -0 elements to 0