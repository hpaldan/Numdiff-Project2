%test_schrodinger

%Different potential equations for testing.
V_1 = @(x) x*0;
V_2 = @(x) 700*(0.5 - abs(x-0.5));
V_3 = @(x) 800*sin(x*pi).^2; 
N = 100;

[eig_modes, eig_vals, dens] = Schrodinger(V_3,N);

%% Plotting eigenfunctions and their corresponding energy level.
x = linspace(0,1,N+2);
for i = 1:1:N
eig_vals(i);
%We want to plot the ones of smallest absolut valie.
plot(x,eig_modes(:,i));
hold on

end

%% Plotting probability densities:
x = linspace(0,1,N+2);
for i = 1:1:N
plot(x,dens(:,i));
hold on

end