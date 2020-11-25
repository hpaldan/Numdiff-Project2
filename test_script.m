
%Numdiff projekt 2 test script
clear all 
close all
clc

N = 100;
L = 4*pi;
alpha = 3;
beta = -7;

f = @(x) cos(x);
x = linspace(0, L, N);
fvec = f(x);

y = twopBVP(fvec, alpha, beta, L, N);
y(1) = alpha;
y(end) = beta;
y
plot(x,y);