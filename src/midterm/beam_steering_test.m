%% Beam steering test
%% script for mid-term
close all
clc
clear
%% Parameters
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;
N = 1000;                % to have a aperture of 1020 lambda
d = param.lambda;
L = 1000;             % now Lz = 700 or 900 doesn't work
win.leftend = -18;
win.rightend = 18;
Res = (win.rightend - win.leftend)/(L - 1);
angle = linspace(win.leftend,win.rightend,L)';  % in degree
angle_90 = linspace(-90,90,5*L)';
theta = angle*pi/180;
theta_90 = angle_90*pi/180;
u = linspace(sin(win.leftend*pi/180),sin(win.rightend*pi/180),L)';
x = (-(N-1)*d/2:d:(N-1)*d/2)';
% x = (0:d:(N-1)*d)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
S_u = exp(1i*param.k*u*x');
S_90 = exp(1i*param.k*sin(theta_90)*x');            
%% Set the desired beam pattern (Loose than our requirement) 
% left_win = -0.4;
% right_win = 0.4;
left_win = -0.04;
right_win = 0.04;

desired_pattern = upperboundgen(win.leftend,win.rightend,[left_win right_win],L);       % The desired pattern
figure
plot(angle,desired_pattern)
xlabel('Angle (degrees)');
ylabel('Amplitude')
title('Ideal beam pattern')
%% Modify the S matrix and reference pattern
loose_range = 0.02;    % in degree   (0.014 is not feasible)
% loose_range = 0.2;

desired_pattern_modified = adjustrefpattern(desired_pattern,[win.leftend win.rightend],[left_win right_win],L,loose_range);

S_new = adjustSmatrix(S,[win.leftend win.rightend],[left_win right_win],L,loose_range);
%% Instead of runnning the following sections, load the result directly. The result is stored at the "W_4_midterm.mat"
load("W_4_midterm.mat")

%% 
figure
subplot(2,1,1)
plot(angle,20*log10(abs(S*w)))
subplot(2,1,2)
plot(angle,20*log10(abs(S*W_rew(:,5))))

%% steering beamformer
angle_t = 10;
theta_t = angle_t*pi/180;
w_steering = exp(-1i*param.k*x*sin(theta_t));

%%
figure
subplot(2,1,1)
w1 = w.*w_steering;
plot(angle,20*log10(abs(S*w1)))
subplot(2,1,2)
w2 = W_rew(:,5).*w_steering;
plot(angle,20*log10(abs(S*w2)))
%%
figure
w_all1 = ones(N,1);
w_matched = exp(-1i*param.k*x*sin(5));      % actually all 1 for broadside direction
subplot(2,1,1)
% plot(angle,20*log10(abs(S*w_all1)))
plot(angle,abs(S*w_all1))
subplot(2,1,2)
% plot(angle,20*log10(abs(S*w_matched)))
plot(angle,abs(S*w_matched))