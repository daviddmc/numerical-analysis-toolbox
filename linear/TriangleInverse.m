function [ A ] = TriangleInverse( A , shape)
% TriangleInverse     Inverse of triangular matrix.
%   TriangleInverse(A) is the inverse of the square lower triangular matrix 
%   A
%
%   TriangleInverse(A, 'upper') is the inverse of the square upper
%   triangular matrix A.
%
%   See also Inverse

%   Copyright 2017 Junshen Xu

CheckSquareMatrix(A, 'A');
n = size(A,1);

if(~exist('shape', 'var'))
    shape = 'lower';
end

A(1:n+1:end) = 1./A(1:n+1:end);

if strncmpi(shape, 'u', 1)
    %{
    for i = 1 : n
        for j = i+1:n
            A(i,j) = A(i,i:j-1) * A(i:j-1, j) * -A(j,j);
        end
    end
    %}
    for j = n : -1 : 2
        for i = j-1: -1 :1
            A(i,j) = A(i,i+1:j) * A(i+1:j, j) * -A(i,i);
        end
    end
elseif strncmpi(shape, 'l', 1)
    for i = n : -1 : 2
        for j = i-1: -1 :1
            A(i,j) = A(i,j+1:i) * A(j+1:i, j) * -A(j,j);
        end
    end
else
    error('shape should be ''upper'' or ''lower''.');
end

