function [ output_args ] = RosenbrockSearch( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


n = length(x0);
A = eye(n);
for iter = 1 : maxIter
    for i = 1:n
        e = OneStep(fun, A(:,i));
        A(:,i) = A(:,i) * e;
    end
    A = cumsum(A, 2, 'reverse');
    x = x + A(:, 1);
    f = fun(x);
    [A, ~] = QR(A);
end


end

function e = OneStep(fun, d, alpha, beta, x0, e0)
e = e0;
f0 = fun(x0);
f = fun(x0 + e * d);
if f < f0
    e = e * alpha;
else
    e = e * beta;
end

end
