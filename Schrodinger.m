%Schrödinger
function [eig_modes, eig_vals, dens] = Schrodinger(V,N)

%We need to create distance dx.
dx = 1/(N+1);

%We need to divide the potential function into a vector. Evaluated at x_i
%where i goes from 1->N. and x_i from 0+dx -> 1-dx. Then the boundary
%values will settle in and decide psi(0) and psi(1). 
%We need to create the modified Toeplitz matrix.
%The modified Toeplitz matrix can consist of 2 parts. One part that is the
%recognized Toeplitz matrix. The other is based on the potential V,
%since we want to go from Psi''-V*psi = -E*psi to T_dx*psi = -E*psi.

%Toeplitz matrix (Same as previous assignment):
R = [-2 1 zeros(1,N-2)];
toep = toeplitz(R);

%Second matrix derived from -V*psi:
V_vec = [zeros(1,N)];
for i = 1:N
    V_vec(i) = V(dx*i);
end
V_mat = diag(V_vec);
V_mat = V_mat;

%Adding the two matrices together and scale them with 1/dx^2
% Jag är lite osäker på om de ska skalas med 1/dx^2 här, värt att
% dubbelkolla en extra gång.
T = 1/dx^2*(toep-V_mat); 

[modes, eig_temp] =  eig(T);

%Puts all the eigenvalues on a vector instead of a matrix.
eigs = diag(eig_temp);
eigs = eigs';
[eigs ind] = sort(eigs,'descend');
modes = modes(:,ind);

%Adding the boundary conditions.
first_line = zeros(1,N);
end_line = zeros(1,N);
modes = [first_line; modes; end_line];
clear first_line end_line

%Rewriting the eigenvalues and eigenmodes and multiplying with -1 to
%accommodate for the -1 in -E*psi.
eig_vals= -1*eigs;
eig_modes = -1*modes;

%Getting the probability density:
dens = zeros(N+2,N);
for i = 1:N
dens(:,i) = abs(eig_modes(:,i)).^2;
end
clear eigs modes








