%% The Beam Equation
clear all
close all 
clc


alpha = 0; %boundary value at x=0
beta = 0; %boundary value at x=L

N = 100;

L = 10; %length of beam
E = 1.9*10^(11); %construction material is steel

x = linspace(0, L, N);


I = @(x) 10.^(-3)*(3-2*(cos(x*pi/L)).^12);
Ivec = I(x);

q = -50000; %load on the beam
qvec = q*ones(size(x));

M = twopBVP(qvec, alpha, beta, L, N);
M = M';

Bendf = M./(E*Ivec);

u = twopBVP(Bendf, alpha, beta, L, N);

%plot(x,M)
plot(x,u)
min(u)

