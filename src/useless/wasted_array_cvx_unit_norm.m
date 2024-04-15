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
desired_pattern = upperboundgen(leftend,rightend,[left_win right_win],L);       % The desired pattern
figure
plot(angle,desired_pattern)
%%
% nc = 500;
    cvx_begin
        variable w(N) complex
        variable x(N) complex   
        minimize(square(norm(x.*w - 1 ,2)) + lambda*norm(w,1))
        subject to
            % 20*log10(abs(S*w/nnz(w)))               % N is the # of active elements
            % 20*log(abs(S*w/nnz(w))/log(10))
            % not good, having problem here.......
            norm(S*w - desired_pattern,2) <= epsilon_;       % problem, not in dB
            w = conj(x);
    cvx_end
%% Calculate and plot the l1 norm of the complex weights
l1n = zeros(nc);            % row vector
for i = 1:1:nc
    l1n(i) = norm(W(:,i),0);
end
epsilon = 1:1:nc;           % row vector
figure;
plot(epsilon,l1n)
xlabel('epsilon')
ylabel('l1 norm')

%% 
% nc = 500, a plot every 20. # of plot = 25
% plot # = 20*i
figure
for i = 1:1:25
    subplot(5,5,i);
    plot(abs(W(:,20*i)))
    title(['No.', num2str(20*i)])
end
sgtitle('Amplitude distribution of complex weights of No.20, 40,..., 500')

%%
figure
for i = 1:1:25
    subplot(5,5,i)
    plot(abs(W(:,2*i)))
    title(['No.', num2str(2*i)])
end
sgtitle('Amplitude distribution of complex weights of No.2, 4,..., 50')

%%
figure
for i = 1:1:25
    subplot(5,5,i)
    plot(angle,20*log10(abs(S*W(:,20*i))))
    title(['No.', num2str(20*i)])
end
sgtitle('Array pattern with complex weights of No.20, 40,..., 500')
%%
figure
for i = 1:1:25
    subplot(5,5,i)
    plot(angle,20*log10(abs(S*W(:,2*i))))
    title(['No.', num2str(2*i)])
end
sgtitle('Array pattern with complex weights of No.2, 4,..., 50')
%% Modify the S matrix and reference pattern
S_new = adjustSmatrix(S,leftend,rightend,[left_win right_win],L,3);
desired_pattern_modified = adjustrefpattern(desired_pattern,leftend,rightend,[left_win right_win],L,3);
%%
% nc = 500;
% W_modified = zeros(N,nc);
% for i = 1:1:nc
    cvx_begin
        variable w(N) complex
        minimize(norm(w,1))
        subject to
            norm(S_new*w - desired_pattern_modified,2) <= 0.01;       % problem, not in dB
    cvx_end
    % W_modified(:,i) = w;
% end

%% Plot the distribution of the weights' amplitude
plot(angle,20*log10(abs(S*w)))
figure;
plot(abs(w1))
% choosing a value as threshold, counting the # of nonzero element
threshold_0 = 0.05;
nnz(max(abs(w1)-threshold_0,0))
