%% This script exam the choice of epsilon
close all
clc
clear
load('epsilon_investigation.mat')       % load the stored data and directly plot graph without running cvx
%% Parameters
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;
N = 1000;                % to have a aperture of 1020 lambda
d = param.lambda;
L = 1000;             % The resolution is at least 0.035 degree
fov.left = -18;
fov.right = 18;
Res = (fov.right - fov.left)/(L - 1);
angle = linspace(fov.left,fov.right,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
Sn = S./sqrt(N);
%% Set the desired beampattern
beam.left = -0.04;
beam.right = 0.04;
desired_pattern = upperboundgen(fov.left,fov.right,[beam.left beam.right],L);
%% Modify the S matrix and reference pattern
desired_pattern_modified = desired_pattern;
S_modified = S;
% loose_range = 0.005;
% S_modified = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
% desired_pattern_modified = adjustRefpattern2(desired_pattern,[fov.left fov.right],[beam.left beam.right],L,loose_range);
%%
ep_log = [0.1 1 10 100];      % 0.1 is infeasible
% ep_log = [1 10 100];
W_log_increase = zeros(N,size(ep_log,2));
%% investigate the right epsilon
for i = 1:1:size(ep_log,2)
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S_modified*w - desired_pattern_modified,2) <= ep_log(i);
cvx_end
W_log_increase(:,i) = w;
end
%%
figure
for i = 1:1:size(ep_log,2)
    subplot(2,2,i)
    plot(angle,20*log10(abs(S_modified*W_log_increase(:,i))))
    xlabel('Angle in degree')
    ylabel('Intensity in dB')
    title("relaxation of epsilon =" + ep_log(i))
end
% conclusion: the best range of epsilon at this particular setup is around
% 1
%%
ep = 1:-0.2:0.2;
W_decrease = zeros(N,size(ep,2));
for i = 1:1:size(ep,2)
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S_modified*w - desired_pattern_modified,2) <= ep(i);
cvx_end
W_decrease(:,i) = w;
end
%%
figure
for i = 1:1:size(ep,2)
    subplot(2,3,i)
    plot(angle,20*log10(abs(S_modified*W_decrease(:,i))))
    xlabel('Angle in degree')
    ylabel('Intensity in dB')
    title("relaxation of epsilon =" + ep(i))
end
% conclusion: the best range for epsilon at thie particular setup is
% between 0.4 and 0.2
%% 
ep_increase = 0.2:0.05:0.4;
W_increase = zeros(N,size(ep_increase,2));
for i = 2:1:size(ep_increase,2)
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S_modified*w - desired_pattern_modified,2) <= ep_increase(i);
cvx_end
W_increase(:,i) = w;
end
%%
figure
for i = 2:1:size(ep_increase,2)
    subplot(2,2,i-1)  
    plot(angle,20*log10(abs(S_modified*W_increase(:,i))))
    xlabel('Angle in degree')
    ylabel('Intensity in dB')
    title("relaxation of epsilon =" + ep_increase(i))
end
% conclusion: the optimal epsilon is in (0.35,0.4]
%%
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S_modified*w - desired_pattern_modified,2) <= 0.39;
cvx_end
% conclusion: the tightest epsilon is 0.4
%% back to log scale examination, exam the pattern and the amplitude of exciations
figure
for i = 2:1:size(ep_log,2)
    subplot(2,3,i-1)
    plot(angle,20*log10(abs(S_modified*W_log_increase(:,i))))
    ylim([-260 0])
    xlabel('Angle in degree')
    ylabel('Intensity in dB')
    title("relaxation of epsilon =" + ep_log(i))
    subplot(2,3,i+2)
    stem(abs(W_log_increase(:,i)))
    ylim([0 6e-4])
    xlabel('Index')
    ylabel('Amplitude')
    title("relaxation of epsilon =" + ep_log(i))
end
%% Since epsilon = 0.4 is the smallest value for feasible problem
L = 3000;             % The resolution is at least 0.035 degree
fov.left = -18;
fov.right = 18;
Res = (fov.right - fov.left)/(L - 1);
angle = linspace(fov.left,fov.right,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
figure
plot(angle,20*log10(abs(S*W_increase(:,end))))
title("The beampattern with sparsest solution and L="+L)