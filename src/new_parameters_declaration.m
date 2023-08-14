%% Array thinning cvx, unfinished
close all
clc
clear
%% Parameters
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;
N = 100;                % to have a aperture of 1020 lambda
d = 0.5*param.lambda;
L = 3000;             % The resolution is at least 0.035 degree
win.leftend = -90;
win.rightend = 90;
Res = (win.rightend - win.leftend)/(L - 1);
angle = linspace(win.leftend,win.rightend,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix

%% Set the desired beam pattern
left_win = -10;
right_win = 10;
desired_pattern = upperboundgen(win.leftend,win.rightend,[left_win right_win],L);       % The desired pattern
figure
plot(angle,desired_pattern)

%% Modify the S matrix and reference pattern
S_new = adjustSmatrix(S,win.leftend,win.rightend,[left_win right_win],L,0.005);
desired_pattern_modified = adjustrefpattern(desired_pattern,win.leftend,win.rightend,[left_win right_win],L,0.005);
