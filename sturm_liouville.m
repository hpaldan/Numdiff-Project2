%Sturm_Liouville 

%Defining all constants
N = 400;
alph = 0;
bet = 1;
dx = 1/(N+1);

%Creating the grid
f_vec = linspace(alph,bet,N+1);

%From lecture slide 22/67
Y_last = 1/3(2*bet*dx+4*Y_f-Y_ff);

R = [-2 1 zeros(1,N-2)];
toep = 1/dx^2* toeplitz(R);
[mode, A] =  eig(toep);

eig_vals = [zeros(1,length(A))];

string_vec = [zeros(1,length(A))];



%%
for i = 1:length(A)

    plot(mode(:,i));
    hold on

end



