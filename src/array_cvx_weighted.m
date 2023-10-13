%% This script try to investigate the correct normalization for S matrix
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
u = linspace(sin(win.leftend*pi/180),sin(win.rightend*pi/180),L)';
x = (-(N-1)*d/2:d:(N-1)*d/2)';
% x = (0:d:(N-1)*d)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
S_u = exp(1i*param.k*u*x');
S_normalized = S./sqrt(N);

%% Set the desired beam pattern (Loose than our requirement) 
% left_win = -0.4;
% right_win = 0.4;
left_win = -0.04;
right_win = 0.04;

desired_pattern = upperboundgen(win.leftend,win.rightend,[left_win right_win],L);       % The desired pattern
figure
plot(angle,desired_pattern)
nzindex = find(desired_pattern == 0);

desired_pattern_with_noise = desired_pattern + wgn(L,1,-48);
desired_pattern_with_margin_for_noise = desired_pattern;
desired_pattern_with_margin_for_noise(nzindex) = desired_pattern_with_margin_for_noise(nzindex) + 2e-2;
%% Modify the S matrix and reference pattern
loose_range = 0.02;    % in degree   (0.014 is not feasible)
% loose_range = 0.2;
S_newn = adjustSmatrix(S_normalized,[win.leftend win.rightend],[left_win right_win],L,loose_range);
desired_pattern_modified = adjustrefpattern(desired_pattern,[win.leftend win.rightend],[left_win right_win],L,loose_range);

desired_pattern_with_noise_modified = adjustrefpattern(desired_pattern_with_margin_for_noise,[win.leftend win.rightend],[left_win right_win],L,loose_range);

S_new = adjustSmatrix(S,[win.leftend win.rightend],[left_win right_win],L,loose_range);
%% cvx gives result with 1/sqrt(N)
cvx_begin
    variable wn(N) complex
    minimize(norm(wn,1))
    subject to
        norm(S_newn*wn - desired_pattern_modified,2) <= 0.1;       % 0.1 for L=1000, 0.72 for L=10000
        % norm(S_newn*wn - desired_pattern_with_noise_modified,2) <= 0.1;      
cvx_end
%% cvx 
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S_new*w - desired_pattern_modified,2) <= 0.1;       % 0.1 for L=1000, 0.72 for L=10000
        % norm(S_newn*wn - desired_pattern_with_noise_modified,2) <= 0.1;      
cvx_end

%% weighted l1 norm
Z = eye(N);
P = 5;
W_rew = zeros(N,P);
for i = 1:1:P    
cvx_begin
    variable w_rew(N) complex
    minimize(norm(Z*w_rew,1))
    subject to
        norm(S_new*w_rew - desired_pattern_modified,2) <= 0.2;       % 0.1 for L=1000, 0.72 for L=10000
        % norm(S_newn*wn - desired_pattern_with_noise_modified,2) <= 0.1;      
cvx_end
W_rew(:,i) = w_rew;
Z = inv(diag(abs(w_rew) + 0.0002));
end
%% plot the weighted result

figure;
for i = 1:1:P
subplot(P,1,i)
semilogy(flip(sort(abs(W_rew(:,i)))))
end
%%
figure;
for i = 1:1:P
subplot(P,1,i)
plot(angle,20*log10(abs(S*W_rew(:,i))))
end

%%
figure;
plot(angle,20*log10(abs(S_normalized*W_rew(:,1))))
hold on
plot(angle,20*log10(abs(S_normalized*W_rew(:,2))))
plot(angle,20*log10(abs(S_normalized*W_rew(:,5))))
% subplot(3,1,2)
% semilogy(flip(sort(abs(W_rew(:,2)))))
% 
% subplot(3,1,3)
% semilogy(flip(sort(abs(W_rew(:,3)))))
%%
L = 50000;
win.leftend = -90;
win.rightend = 90;
angle = linspace(win.leftend,win.rightend,L)';  % in degree
theta = angle*pi/180;
u = linspace(sin(win.leftend*pi/180),sin(win.rightend*pi/180),L)';
x = (-(N-1)*d/2:d:(N-1)*d/2)';
% x = (0:d:(N-1)*d)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
S_u = exp(1i*param.k*u*x');
S_normalized = S./sqrt(N);
figure;
plot(angle,20*log10(abs(S_normalized*W_rew(:,5))))

%% cvx with reference pattern left up by noise
cvx_begin
    variable wn(N) complex
    minimize(norm(wn,1))
    subject to
        % 0.1 for L=1000, 0.72 for L=10000
        norm(S_newn*wn - desired_pattern_with_noise_modified,2) <= 0.1;      
cvx_end

%% ISTA
S_tilda = [real(S) -imag(S);imag(S) real(S)];
desired_pattern_tilda = [real(desired_pattern);imag(desired_pattern)];
%%
[w_ista,err] = ISTA(desired_pattern_tilda,S_tilda,1,100,N);
w_ISTA = complex(w_ista(1:N),w_ista(N+1:end));
figure
subplot(3,1,1)
plot(err)
subplot(3,1,2)
stem(abs(w_ISTA))
subplot(3,1,3)
plot(angle,20*log10(abs(S*w_ISTA)))

%% OMP
S_tilda_OMP = [real(S_normalized) -imag(S_normalized);imag(S_normalized) real(S_normalized)];
w_omp = OMP(desired_pattern_tilda,S_tilda_OMP,0.1,N);
w_OMP = complex(w_omp(1:N),w_omp(N+1:end));
figure
subplot(2,1,1)
stem(abs(w_OMP))
subplot(2,1,2)
plot(angle,20*log10(abs(S*w_OMP)))

%%
disp("with normalized S, cvx gives a result w with l2 norm of " +norm(wn))

%%
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S_new*w - desired_pattern_modified,2) <= 0.1;       % (problem, not in dB) not really a problem
cvx_end
%%
% [cost,th_level,w_th,W,MSE,MSE_new,MSE_app] = iterativethreshold(wn,S,0.002);
[cost,th_level,w_th,W,MSE,MSE_new,MSE_app] = iterativethreshold(w,S,0.002);

%% Result evaluation (x axis is sparsity)
k = 0:1:N-1;

nnze = zeros(N,1);
for i = 1:1:N
    % MSE(i) = (W(:,i) - w)'*S'*S*(W(:,i) - w)/L;
    nnze(i) = nnz(W(:,i));
end
subplot(3,1,1)
semilogy(k,flip(abs(MSE)))
% plot(k,flip(abs(MSE)))
xlabel('Sparsity')
ylabel('MSE')
title('The MSE of the threshold pattern vs sparsity')

subplot(3,1,2)
semilogy(k,flip(abs(MSE_app)))
% plot(k,flip(abs(MSE_app)))
xlabel('Sparsity')
ylabel('MSE approximation')
title('The MSE of the threshold pattern vs sparsity')

subplot(3,1,3)
plot(k,flip(nnze))
xlabel('Sparsity')
ylabel('Number of non-zeros entries')
title('The number of nonzero element vs sparsity')

%%

function [z,err] = ISTA(y,A,lambda,iter,N)
eigenlambda = 2*max(eig(A'*A));  % Lippschitz can be anything larger than bigest eigenvalue
eta = 1/eigenlambda;         % inverse of that
shrink_threshold = lambda*eta;         % the threshold 
z = ones(2*N,1);
err = zeros(iter,1);
% iteration
for i = 1:1:iter 
    p = z - 2*eta*A'*(A*z - y);     % gradient step
    z = (max(abs(p) - shrink_threshold,0)).*sign(p);
    err(i) = mean(norm(A*z - y,2)); 
end
end

%% Function OMP
% Given "measurement" y, "CS matrix" A, "sparsity" k
% function [z,err] = OMP(y,A,epsilon)
function z = OMP(y,A,epsilon,N)
support=[];
r=y;  
z = zeros(2*N,1);
% prod = zeros(N,1);
for times = 1:1:10000
    prod = abs(A'*r);
    [val,pos] = max(prod);
    support = [support, pos];
    % LS
    z_lambda = pinv(A(:,support))*y;
    r = y - A(:,support)*z_lambda;
%     err_(times) = norm(r,2);
    err_ = norm(r,2)^2;
    if err_<epsilon*10 %
        err_
            break; 
    end
end
% err = err_';
z(support) = z_lambda;
end

% %%
% cvx_begin
%     variable w2(N) complex
%     minimize(norm(w2,1))
%     subject to
%         norm(S_new*w2 - desired_pattern_modified,2) <= 2;       % (problem, not in dB) not really a problem
%         w2'*w2 <= 1;
% cvx_end
% disp("with w'*w <= 1, cvx gives a result w with l2 norm of " +norm(w2))
% 
% %%
% disp("Epsilon of 2 gives a meaningful result here, with l2 norm of w as " + norm(w2) +". (3 doesn't work)")
% 
% stem(abs(w2))
% title('Not sparse any more')