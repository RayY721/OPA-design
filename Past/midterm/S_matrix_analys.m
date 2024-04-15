%% script for mid-term
close all
clc
clear
%% Parameters
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;
d = param.lambda;

N = 300;
L = 400;
win.leftend = -10;
win.rightend = 10;

Res = (win.rightend - win.leftend)/(L - 1);
angle = linspace(win.leftend,win.rightend,L)';  % in degree
theta = angle*pi/180;
u = linspace(sin(win.leftend*pi/180),sin(win.rightend*pi/180),L)';
x = (-(N-1)*d/2:d:(N-1)*d/2)';

% x = (0:d:(N-1)*d)';

% x = (0:d:(N-1)*d)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
S_u = exp(1i*param.k*u*x');
%% Display the real part of S matrix
[X,Y] = meshgrid(1:N,1:L);
figure
subplot(2,2,1)
surf(X,Y,real(S_u));
title('The real part of S_u matrix with L numebr of samples over u')
subplot(2,2,2)
surf(X,Y,real(S));
title('The real part of S matrix with L numebr of samples over theta')
DFT_matrix = dftmtx(N);
% IDFT_matrix = conj(dftmtx(N))/N;
IDFT_matrix = conj(dftmtx(N));
[X,Y] = meshgrid(1:N,1:N);
subplot(2,2,3)
surf(X,Y,real(DFT_matrix));
title('The real part of DFT matrix')
subplot(2,2,4)
surf(X,Y,real(IDFT_matrix));
title('The real part of IDFT matrix')
%% Display the S^H * S 
[X,Y] = meshgrid(1:N,1:N);
DFT_matrix = dftmtx(N);
% IDFT_matrix = conj(dftmtx(N))/N;
IDFT_matrix = conj(dftmtx(N));
figure
subplot(2,2,1)
surf(X,Y,real(S_u'*S_u))
subplot(2,2,2)
surf(X,Y,real(S'*S))
subplot(2,2,3)
surf(X,Y,real(DFT_matrix'*DFT_matrix))
subplot(2,2,4)
surf(X,Y,real(IDFT_matrix'*IDFT_matrix))
%%
figure
subplot(2,1,1)
plot(imag(exp(1i*sin(theta))))
subplot(2,1,2)
plot(imag(exp(1i*u)))
% u and theta have barely difference as the FOV is confined 
% to [-18,18], where the 
