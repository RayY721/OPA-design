close all
clc
clear

%% Parameters
% useless parameters
param.c = physconst('LightSpeed');

% useful parameters
param.lambda = 1550e-9;
param.k = 2*pi/param.lambda;
N = 1000;
spacing = param.lambda;
L = 1000;
fov.leftend = -18;
fov.rightend = 18; 
angle = linspace(fov.leftend,fov.rightend,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*spacing/2:spacing:(N-1)*spacing/2)';
% S and normalized S matrix
S = exp(1i*param.k*sin(theta)*x');            % S matrix
Sn = S./sqrt(N);
%% Set the desired beam pattern
left_win = -0.04;
right_win = 0.04;   % overall beamwidth is 0.08 degree
desired_pattern = upperboundgen(fov.leftend,fov.rightend,[left_win right_win],L);
figure
plot(angle,desired_pattern)
xlabel('Angle (degrees)');
ylabel('Amplitude')
title('Ideal beam pattern')
%% Modify the S matrix and reference pattern
loose_range = 0.02;    % in degree   (0.014 is not feasible)
% loose_range = 0.2;
S_new = adjustSmatrix(S,[fov.leftend fov.rightend],[left_win right_win],L,loose_range);
Sn_new = adjustSmatrix(Sn,[fov.leftend fov.rightend],[left_win right_win],L,loose_range);
desired_pattern_modified = adjustrefpattern(desired_pattern,[fov.leftend fov.rightend],[left_win right_win],L,loose_range);

%% CVX with S matrix epsilon 0.1    % not converging, no solution
cvx_begin
    variable w1(N) complex
    minimize(norm(w1,1))
    subject to
        norm(S*w1 - desired_pattern,2) <= 0.1;
cvx_end
%% CVX with Sn matrix epsilon 0.1   % not converging, no solution
cvx_begin
    variable w2(N) complex
    minimize(norm(w2,1))
    subject to
        norm(Sn*w2 - desired_pattern,2) <= 0.1;
cvx_end
%% CVX with S_new matrix epsilon 0.1    % not converging, no solution
cvx_begin
    variable w3(N) complex
    minimize(norm(w3,1))
    subject to
        norm(S_new*w3 - desired_pattern_modified,2) <= 0.1;
cvx_end
%% CVX with Sn_new matrix epsilon 0.1   % not converging, no solution
cvx_begin
    variable w4(N) complex
    minimize(norm(w4,1))
    subject to
        norm(Sn_new*w4 - desired_pattern_modified,2) <= 0.1;
cvx_end
%%
figure
subplot(2,2,1)
plot(angle,20*log10(abs(S*w1)))
subplot(2,2,2)
plot(angle,20*log10(abs(Sn*w2)))
subplot(2,2,3)
plot(angle,20*log10(abs(S*w3)))
subplot(2,2,4)
plot(angle,20*log10(abs(Sn*w4)))
%% CVX with S matrix epsilon 0.2    % 
cvx_begin
    variable w5(N) complex
    minimize(norm(w5,1))
    subject to
        norm(S*w5 - desired_pattern,2) <= 0.1;
cvx_end
%% CVX with Sn matrix epsilon 0.2   % 
cvx_begin
    variable w6(N) complex
    minimize(norm(w6,1))
    subject to
        norm(Sn*w6 - desired_pattern,2) <= 0.1;
cvx_end
%% CVX with S_new matrix epsilon 0.2    % 
cvx_begin
    variable w7(N) complex
    minimize(norm(w7,1))
    subject to
        norm(S_new*w7 - desired_pattern_modified,2) <= 0.1;
cvx_end
%% CVX with Sn_new matrix epsilon 0.2   % 
cvx_begin
    variable w8(N) complex
    minimize(norm(w8,1))
    subject to
        norm(Sn_new*w8 - desired_pattern_modified,2) <= 0.1;
cvx_end