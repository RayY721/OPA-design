%% This script try to investigate relaxation allocation over angular domain
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
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
S_normalized = S./sqrt(N);

%% Set the desired beam pattern (Loose than our requirement) 
% left_win = -0.4;
% right_win = 0.4;
left_win = -0.04;
right_win = 0.04;

desired_pattern = upperboundgen(win.leftend,win.rightend,[left_win right_win],L);       % The desired pattern
figure
plot(angle,desired_pattern)
% desired_pattern_with_noise = desired_pattern + wgn(L,1,-48);
%% Modify the S matrix and reference pattern
loose_range = 0.02;    % in degree   (0.014 is not feasible)
% loose_range = 0.2;
S_new = adjustSmatrix(S_normalized,[win.leftend win.rightend],[left_win right_win],L,loose_range);
desired_pattern_modified = adjustrefpattern(desired_pattern,[win.leftend win.rightend],[left_win right_win],L,loose_range);

% S_new = adjustSmatrix(S,[win.leftend win.rightend],[left_win right_win],L,loose_range);
%%
w_vector = zeros(N,10);
ep = 0.1:0.1:1;
for i = 1:1:10
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S_new*w - desired_pattern_modified,2) <= ep(i);       % 0.1 for L=1000, 0.72 for L=10000
cvx_end
w_vector(:,i) = w;
end

%%
for i = 1:1:10
    % norm(w_vector(:,i)) 
    w_vector(:,i)'*w_vector(:,i)
end

%%
figure
for i = 1:1:10
    subplot(5,2,i)
    plot(angle,20*log10(abs(S_normalized*w_vector(:,i))))
    xlabel('Angle in degree')
    ylabel('Intensity in dB')
    title("relaxation of epsilon =" + ep(i))
end
%%
figure
for i = 1:1:10
    subplot(5,2,i)
    stem(abs(desired_pattern_modified - abs(S_new*w_vector(:,i))))
    ylim([0 0.8])
end

%%
figure
for i = 1:1:10
    subplot(5,2,i)
    stem(abs(desired_pattern - abs(S_normalized*w_vector(:,i))))
    ylim([0 0.8])
end
%%
w_vector_log_increase = zeros(N,4);
ep_log = [0.1 1 10 100];
for i = 1:1:4
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S_new*w - desired_pattern_modified,2) <= ep_log(i);       % 0.1 for L=1000, 0.72 for L=10000
cvx_end
w_vector_log_increase(:,i) = w;
end

%%
figure
for i = 1:1:4
    subplot(2,2,i)
    plot(angle,20*log10(abs(S_normalized*w_vector_log_increase(:,i))))
    xlabel('Angle in degree')
    ylabel('Intensity in dB')
    title("relaxation of epsilon =" + ep_log(i))
end

%% cvx Lasso
cvx_begin
    variable w2(N) complex
    minimize(norm(S_new*w2 - desired_pattern_modified,2) + 0.01*norm(w2,1))
cvx_end