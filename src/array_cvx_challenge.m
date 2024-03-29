%% Array thinning cvx, challenge
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
L = 700;             % The resolution is at least 0.035 degree
win.leftend = -18;
win.rightend = 18;
Res = (win.rightend - win.leftend)/(L - 1);
angle = linspace(win.leftend,win.rightend,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix

%% Set the desired beam pattern
left_win = -0.04;
right_win = 0.04;
desired_pattern = upperboundgen(win.leftend,win.rightend,[left_win right_win],L);       % The desired pattern
figure
plot(angle,desired_pattern)
%% Modify the S matrix and reference pattern
loose_range = 0.005;    % in degree
S_new = adjustSmatrix(S,[win.leftend win.rightend],[left_win right_win],L,loose_range);
desired_pattern_modified = adjustrefpattern(desired_pattern,[win.leftend win.rightend],[left_win right_win],L,loose_range);

%%
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S_new*w - desired_pattern_modified,2) <= 0.21;       % problem, not in dB
cvx_end
% W_modified_small(:,i) = w;

%% View the pattern over FOV
figure;
subplot(2,2,1)
plot(angle,20*log10(abs(S*w)/max(abs(S*w))))
hold on
xline(9)
xline(-9)
xlabel('Angle in degree')
ylabel('Normalized  intensity')
title('Beam pattern over the optimization range')

% View the whole [-90,90]
Re = 30000;
angle_90 = linspace(-90,90,Re)';  % in degree
theta_90 = angle_90*pi/180;
S_90 = exp(1i*param.k*sin(theta_90)*x');            % S matrix
% figure;
subplot(2,2,2)
plot(angle_90,20*log10(abs(S_90*w)/max(abs(S_90*w))))
hold on
xline(9)
xline(-9)
xlabel('Angle in degree')
ylabel('Normalized  intensity')
title('Beam pattern over the [-90,90] range')

% Pattern steering
w_steering = exp(-1i*param.k*x*sin(5*pi/180));
w_steered = w_steering.*w;
% Using w_steered as an additional beamformer
% figure;
subplot(2,2,3)
plot(angle,20*log10(abs(S*w_steered)/max(abs(S*w_steered))))
hold on
xline(9)
xline(-9)
xlabel('Angle in degree')
ylabel('Normalized  intensity')
title('Steered beam pattern over optimization range')

% figure;
subplot(2,2,4)
plot(angle_90,20*log10(abs(S_90*w_steered)/max(abs(S_90*w_steered))))
hold on
xline(9)
xline(-9)
xlabel('Angle in degree')
ylabel('Normalized  intensity')
title('Steered pattern over [-90,90]')


%% 
Re = 30000;
angle_18 = linspace(-18,18,Re)';  % in degree
theta_18 = angle_18*pi/180;
S_18 = exp(1i*param.k*sin(theta_18)*x');            % S matrix
figure;
plot(angle_18,20*log10(abs(S_18*w_steered)/max(abs(S_18*w_steered))))
hold on
xline(9)
xline(-9)
title('Steered pattern over FOV')

%% Iterative hard thresholding pursuit?

[cost,th_level,w_th,W,MSE,MSE2] = iterativethreshold(w,S,0.002);

%% Result evaluate (x axis is thresholding level)

figure; 
subplot(3,1,1)
plot(sort(abs(w),'ascend'),abs(cost))
hold on
xline(th_level)
title('The cost function vs thresholding level')

% MSE error
% original pattern is the mean(reference), the MSE is the sum of error's
% square and normalized by a factor 1/n.
% MSE = (w_th - w)'*S'*S*(w_th - w)/L

% MSE = zeros(N,1);
nnze = zeros(N,1);
for i = 1:1:N
    % MSE(i) = (W(:,i) - w)'*S'*S*(W(:,i) - w)/L;
    nnze(i) = nnz(W(:,i));
end
subplot(3,1,2)
plot(sort(abs(w),'ascend'),abs(MSE))
hold on
xline(th_level)
title('The MSE of the threshold pattern vs thresholding level')
subplot(3,1,3)
plot(sort(abs(w),'ascend'),nnze)
hold on
xline(th_level)
title('The number of nonzero element vs thresholding level')

%% Result evaluation (x axis is sparsity)
k = 0:1:N-1;
figure; 
subplot(3,1,1)
plot(k,flip(abs(cost)))
xlabel('Sparsity')
ylabel('Value of the objective function')
title('The cost function vs sparsity')

% MSE error
% original pattern is the mean(reference), the MSE is the sum of error's
% square and normalized by a factor 1/n.
% MSE = (w_th - w)'*S'*S*(w_th - w)/L

% MSE = zeros(N,1);
nnze = zeros(N,1);
for i = 1:1:N
    % MSE(i) = (W(:,i) - w)'*S'*S*(W(:,i) - w)/L;
    nnze(i) = nnz(W(:,i));
end
subplot(3,1,2)
plot(k,flip(abs(MSE)))
xlabel('Sparsity')
ylabel('MSE')
title('The MSE of the threshold pattern vs sparsity')

subplot(3,1,3)
plot(k,flip(nnze))
xlabel('Sparsity')
ylabel('Number of non-zeros entries')
title('The number of nonzero element vs sparsity')


%% Find a good threshold to truncate 
figure;
subplot(2,1,1)
stem(abs(w))
ylabel('amplitude')
xlabel('elements index')
subplot(2,1,2)
% semilogy(sort(abs(w)))
stem(abs(w_th))
ylabel('amplitude')
xlabel('elements index')

% %% Truncate the weights
% w_abs = abs(w);
% w_truncate = w;
% w_truncate(w_abs <= 0.0002) = 0;
% nnz(w_truncate)

% %% Measure the error due to the truncation (over the FOV range)
% norm((abs(S*w) - abs(S*w_truncate)),2)^2/L
% 
% % norm(abs(S(w - w_truncate),2)^2

%% Visually check
L = 7000;             % The resolution is at least 0.035 degree
win.leftend = -15;
win.rightend = 15;
Res = (win.rightend - win.leftend)/(L - 1);
angle = linspace(win.leftend,win.rightend,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
figure;
subplot(2,1,1)
plot(angle,20*log10(abs(S*w)/max(abs(S*w))))
xlabel('Angle in degree')
ylabel('Normalized  intensity')
title('Original pattern')
% figure;
subplot(2,1,2)
plot(angle,20*log10(abs(S*w_th)/max(abs(S*w_th))))
xlabel('Angle in degree')
ylabel('Normalized  intensity')
title('Truncated pattern')

%% Following investigation of the scaling amplitude

w_scaled = w./10;
figure;
subplot(2,2,1)
plot(angle,20*log10(abs(S*w)/max(abs(S*w))))
title('Original pattern (normalized)')
% figure;
subplot(2,2,2)
plot(angle,20*log10(abs(S*w_scaled)/max(abs(S*w_scaled))))
title('Scaled pattern (normalized)')
subplot(2,2,3)
plot(angle,20*log10(abs(S*w)))
title('Original pattern (unormalized)')
% figure;
subplot(2,2,4)
plot(angle,20*log10(abs(S*w_scaled)))
title('Scaled pattern (unormalized)')

%% plot many w
figure
for i =1:1:25
    subplot(5,5,i)
    plot(angle,20*log10(abs(S*W(:,20*i))))
end

% for a particular w_t, the normalizing factor is that 1/(s(theta_max)'*w_t)
% it makes the normalized S matrix as S/abs(s(theta_max)'*w_t)

% K-means clustering

% Thresholding with energy percentage
