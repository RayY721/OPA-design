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
N = 500;
L = 400;
win.leftend = -10;
win.rightend = 10;
Res_u = (sin(win.rightend*pi/180) - sin(win.leftend*pi/180))/(L-1);
% Res = (win.rightend - win.leftend)/(L - 1);

angle = linspace(win.leftend,win.rightend,L)';  % in degree
theta = angle*pi/180;
u_ = linspace(sin(win.leftend*pi/180),sin(win.rightend*pi/180),L)';
x = (-(N-1)*d/2:d:(N-1)*d/2)';
% x = (0:d:(N-1)*d)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
S_u = exp(1i*param.k*u*x');
%% Set the desired beam pattern (Loose than our requirement) 
left_win = -0.4;
right_win = 0.4;
% left_win = -0.04;
% right_win = 0.04;
desired_pattern = upperboundgen(win.leftend,win.rightend,[left_win right_win],L);       % The desired pattern
figure
plot(angle,desired_pattern)
xlabel('Angle (degrees)');
ylabel('Amplitude')
title('Ideal beam pattern')
%% for fun, try out the atomic norm
cvx_begin sdp
    variable u(L) complex;
    variable y(L) complex;
    variable t complex;
    T = toeplitz(u);

    minimize real(trace(T) + t)
    subject to 
        norm(desired_pattern - y) <= 1;
        [t y';y T] >= 0;
cvx_end
%%
[V,D] = eig(T*T');
V = V(:,end-rank(T*T')+1:end);
V_u = V(1:end-1,:);
V_l = V(2:end,:);
eigenz = eig(V_u'*V_l,V_u'*V_u);
x = imag(eigenz)./(param.k*Res_u);

%%
% https://github.com/JAMES-YI/00_Separation-free-Super-resolution/blob/main/SANM_cvx.m


% function [rx,err]=SANM_cvx(ox,K,N,kk,stepsize,epsilon)
% %% Standard Atomic Norm Minimization
% 
% % solve the dual problem
% cvx_quiet true
% cvx_precision default
% 
% cvx_begin sdp
%   variable rx(N) complex;
%   variable u(N) complex;
%   variable s complex;
%   Q = toeplitz(u);
% 
%   minimize real(trace(Q)/N + trace(s))
%   subject to 
%         [s rx'; rx Q] >= 0,
%         rx(K)==ox(K);
% cvx_end
% 
% err = norm(rx-ox)/norm(ox);
% end