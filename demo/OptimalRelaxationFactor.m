
%% parameters
n = 9; % size of grid
epsilon = 1e-6; % tolerance of error

%% setup
maxIter = ceil(-log10(epsilon) * log(10) * (n+1)^2 / pi^2); % max # of iteration
% calculate maxtrix A
block = 4 * eye(n) + diag(-ones(n-1,1),1) + diag(-ones(n-1,1),-1);
blocks = cell(n, 1);
blocks(:) = {block};
A = blkdiag(blocks{:}) - diag(ones(n^2-n,1),n) -diag(ones(n^2-n,1),-n);
% set omega to the optimal value
omega = 2 / (1 + sin(pi / (n+1)));
% calculate BSOR
D = diag(diag(A));
invDL = Inverse(D + omega * tril(A, -1));
BSOR = invDL * ((1-omega)*D - omega * triu(A, 1)); 
% calculate BJ
BJ = eye(n^2) - D \ A;
% calculate BGS
invDL = Inverse(tril(A, 0));
BGS = eye(n^2) - invDL * A;
% initailize error and x
eSOR = zeros(maxIter, 0);
eGS = zeros(maxIter, 0);
eJ = zeros(maxIter, 0);
xSOR = ones(n^2,1);
xJ = ones(n^2,1);
xGS = ones(n^2,1);

%% iteration
for k = 1 : maxIter
    % calculate error
    eSOR(k) = max(abs(xSOR));
    eGS(k) = max(abs(xGS));
    eJ(k) = max(abs(xJ));
    if(eSOR(k) > epsilon || eGS(k) > epsilon || eJ(k) > epsilon) 
        % errors are small enough
        xSOR = BSOR * xSOR;
        xGS = BGS * xGS;
        xJ = BJ * xJ;
    else
        % end iteration
        eSOR = eSOR(1:k);
        eGS = eGS(1:k);
        eJ = eJ(1:k);
        break
    end
end

%% show results
iter = 1 : length(eJ);
semilogy(iter, eJ, iter, eSOR, iter, eGS);
xlabel('\# of iteration', 'Interpreter','latex')
ylabel('$||x-x^*||_\infty$', 'Interpreter','latex');
legend('Jcobi','SOR','Gauss-Seidel','location','best');


