function [L, D] = Cholesky(A)
% Cholesky decomposition
% L * L^T = A or L * D * L^T = A
% input
% A
% output
% L
% D

CheckSquareMatrix(A, 'A');

