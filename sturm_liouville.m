
%Sturm_Liouville 

%Defining all constants
N = 99;
alph = 0;
bet = 0;
dx = 1/(N+1);
yvec = linspace(0,1,N);
yvec(1) = yvec(1)-alph/dx^2;
%In order to get the Neumann condition right we need to rewrite Y(N-1) and
%Y(N).
yvec(N) = 

((1-1/3)*y_{N-1} + (- 2 +4/3)*y_N + 2/3*bet*dx  )/dx^2



%Calculates the FDM:
R = [-2 1 zeros(1,N-2)];
toep = 1/dx^2 * toeplitz(R);

%Adds the initial boundary condition.
toep(1,1) = toep(1,1)+(alph*1/dx^2);

%Adds the final boundary condition 
%(note that it is a boundary condition on the derivative).
T_f = 1/3*(2*bet*dx+4*toep(N,N)-toep(N,N-1));
T_f = 0;
toep(N,N) = toep(N,N)+(T_f*(1/dx^2));


%Creates the eigenvaluefunctions (modes) and values for 
%the toeplitzmatrix.
[modes, eig_temp] =  eig(toep);
eigs = [zeros(1,N)];

%Puts all the eigenvalues on a vector instead of a matrix.
eigs = diag(eig_temp);
eigs = eigs';
[eigs ind] = sort(eigs,'descend');
modes = modes(:,ind);

clear eig_temp;
%% Creates som cool plots of the modes
for i = N:1:N
eigs(i)
%We want to plot the ones of smallest absolut valie.
plot(modes(:,i));
hold on

end
