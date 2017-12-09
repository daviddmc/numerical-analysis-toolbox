% This script shows the Gershgorin circle theorem which may be used to 
% bound the spectrum of a square matrix. Let A be a complex N-by-N matrix. 
% for i = 1:N, let R_i be the sum of the absolute values of the 
% non-diagonal entries in the i-th row (or column). Let D(A_ii, R_i) be the
% closed disc centered at A_ii with radius R_i. Such a disc is called a
% Gershgorin disc. The theorem states that every eigenvalue of A lies 
% within at least one of the Gershgorin discs.

% Copyright 2017 Junshen Xu

%% init
N = 7; % size of matrix;
% generate a complex matrix randomly
A = rand(N) + rand(N)*1i + diag(3*randperm(N)) + diag(3*randperm(N))*1i;

%% calculate eigenvalues
E = Eigen(A);
Ex = real(E);
Ey = imag(E);

%% show results
scatter(Ex, Ey);
title({'the dots are eigenvalues of A', ...
    'the circles are Gershgorin discs'});
xlabel('real');
ylabel('imaginary')
hold on
for i = 1:N
    r1 = sum(abs(A([1:i-1 i+1:N], i)));
    r2 = sum(abs(A(i, [1:i-1 i+1:N])));
    r = min(r1, r2); % radius of the discs
    x = real(A(i, i));
    y = imag(A(i, i));
    rectangle('Position',[x-r,y-r,2*r,2*r],'Curvature',[1,1],'linewidth',1)
end
axis equal
grid on
hold off
