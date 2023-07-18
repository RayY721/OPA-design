%% Array thinning cvx
close all
clc
clear
%% Parameters
c = physconst('LightSpeed');
lambda = 1550e-9;
fc = c/lambda;
k = 2*pi/lambda;
N = 100;                % to have a aperture of 1020 lambda
d = 0.5*lambda;
L = 3000;             % The resolution is at least 0.035 degree
leftend = -90;
rightend = 90;
Res = (rightend - leftend)/(L - 1);
angle = linspace(leftend,rightend,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';

%%
S1 = exp(transpose(1i*k*x*sin(theta')));   % alternative
S2 = exp(1i*k*sin(theta)*x');
% S and S2 are suppose to be the same???
% But why it's not zero when subtrating them
Sd = S1 - S2;
nnz(Sd)
%% Conclusion
% Two expressions above are the same, the differnece is due to the machine
% precision