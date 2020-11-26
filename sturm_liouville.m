

%Sturm_Liouville 

%Defining all constants
N = 499;
alph = 0;
bet = 0;
dx = 1/(N+1);

%In order to get the Neumann condition right we need to rewrite Y(N-1) and
%Y(N).

%((1-1/3)*y_{N-1} + (- 2 +4/3)*y_N + 2/3*bet*dx  )/dx^2



%Calculates the FDM:
R = [-2 1 zeros(1,N-2)];
toep = 1/dx^2 * toeplitz(R);

%Adds the initial boundary condition.
toep(1,1) = toep(1,1)+(alph*1/dx^2);



%Creates the eigenvaluefunctions (modes) and values for 
%the toeplitzmatrix.
[modes, eig_temp] =  eig(toep);


%Puts all the eigenvalues on a vector instead of a matrix.
eigs = diag(eig_temp);
eigs = eigs';
[eigs ind] = sort(eigs,'descend');
modes = modes(:,ind);

%Adding the final boundary conditions.
end_line = ones(1,N)*(1/3*2*bet*dx)+(1/3*4.*modes(end,:))-(1/3*modes(end-1,:));
first_line = alph*ones(1,N);
modes = [first_line; modes; end_line];

clear eig_temp;
%% Creates som cool plots of the modes
for i = N-5:1:N
eigs(i)
%We want to plot the ones of smallest absolut valie.
plot(modes(:,i));
hold on

end
