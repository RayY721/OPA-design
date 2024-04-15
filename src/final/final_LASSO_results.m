%% Detailed investigation on the evolution of the epsilon. 
%% This script exam the choice of the epsilon
close all
clc
clear
% load('LASSO_results_and_epsilons_used.mat')
    
%% Parameters
for j = 1:1
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;
N = 1100;                % to have a aperture of 1020 lambda
d = param.lambda;
L = 1100;             % The resolution is at least 0.035 degree
fov.left = -18;
fov.right = 18;
Res = (fov.right - fov.left)/(L - 1);
angle = linspace(fov.left,fov.right,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
steering_vector_0 = exp(1i*param.k*sin(0)*x');    
end
%% Set the desired beampattern
beam.left = -0.035;
beam.right = 0.035;
p_d = upperboundgen(fov.left,fov.right,[beam.left beam.right],L);

%%%%%%%%%%%%%%%%%%%%%% LS solution investigation %%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the LS solution with pinv()
w_ls = pinv(S)*p_d;     % This solution gives the LS error
epsilon_ls = norm(S*w_ls - p_d,2);
epsilon_max = norm(p_d,2);

%% CVX for LS solution
cvx_begin
    variable w_ls_cvx(N) complex
    minimize(norm(S*w_ls_cvx - p_d))
cvx_end

epsilon_ls_cvx = norm(S*w_ls_cvx - p_d,2);

%% plot the beampattern and distribution
for j = 1:1
figure
plot(angle,20*log10(abs(S*w_ls)));
set(gca,'FontSize',12);
xlabel('\theta','FontSize',16);
ylabel('Intensity(dB)','FontSize',16);
title('Beampattern','FontSize',16);
end
%% the distribution
for j = 1:1
figure
stem(abs(w_ls))
xlim([0,1100])
set(gca,'FontSize',12);
xlabel('Index of grid points','FontSize',16);
ylabel('Amplitude','FontSize',16);
title('Excitation over the Grid','FontSize',16);
end
%% beampattern
for j = 1:1
figure
plot(angle,20*log10(abs(S*w_ls_cvx)));
set(gca,'FontSize',12);
xlabel('\theta','FontSize',16);
ylabel('Intensity(dB)','FontSize',16);
title('Beampattern','FontSize',16);
end
%% the distribution
for j = 1:1
figure
stem(abs(w_ls_cvx))
xlim([0,1100])
set(gca,'FontSize',12);
xlabel('Index of grid points','FontSize',16);
ylabel('Amplitude','FontSize',16);
title('Excitation over the Grid','FontSize',16);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CVX for varying epsilon
epsilon_1 = linspace(epsilon_ls_cvx,epsilon_max,25);
w_epsilon_increase = zeros(N,size(epsilon_1,2));
% Later on, a non-uniform increase of epsilon is used. 

%% Validate the LASSO problem with epsilon set to LS value
cvx_begin
    variable w_ls_cvx_lasso(N) complex
    minimize(norm(w_ls_cvx_lasso,1))
    subject to
        norm(S*w_ls_cvx_lasso - p_d,2) <= epsilon_ls_cvx;
cvx_end
%% Have a closer look at the 
epsilon_3 = linspace(epsilon_ls_cvx,0.4,25);
w_epsilon_increase3 = zeros(N,size(epsilon_3,2));
%
for i = 1:1:size(epsilon_3,2)
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S*w - p_d,2) <= epsilon_3(i);
cvx_end
w_epsilon_increase3(:,i) = w;
end

%% Plot the spatial distribution of the excitations 
figure
for i =1:M
    subplot(5,5,i)
    stem(abs(w_epsilon_increase3(:,i)))
    xlim([0,1100])
    xlabel('Elements index');
    ylabel('Amplitude');
    title("Excitation distribution \epsilon = " + epsilon_3(i));
end
%% Plot the sorted amplitude distribution of the excitations
figure
for i =1:M
    subplot(5,5,i)
    semilogy(sort(abs(w_epsilon_increase3(:,i)),'descend'))
    xlim([0,1100])
    xlabel('Sorted Index');
    ylabel('Amplitude(dB)');
    title("Excitation Distribution \epsilon = " + epsilon_3(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for i = 1:1:size(epsilon_1,2)
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S*w - p_d,2) <= epsilon_1(i);
cvx_end
w_epsilon_increase(:,i) = w;
end

%% plot the l1 norm of the solutions
[N,M] = size(w_epsilon_increase);
l1_norms = zeros(1,M);
for i = 1:M
    l1_norms(i) = norm(w_epsilon_increase(:,i),1);
end
figure
plot([epsilon_2(2:end) epsilon_1(6:end) ],[l1_norms2(2:end) l1_norms(6:end)],'LineWidth',2);
set(gca,'FontSize',12);
xlabel('Epsilon','FontSize',16);
ylabel('L1 Norm of Excitation Vector','FontSize',16);
title('Norm of Vectors in W Against Epsilon','FontSize',16);
%%
figure
for i =1:M
    subplot(5,5,i)
    stem(abs(w_epsilon_increase(:,i)))
    xlim([0,1100])
    xlabel('Elements index');
    ylabel('Amplitude');
    title("Excitation distribution \epsilon = " + epsilon_1(i));
end
%%
figure
for i =1:M
    subplot(5,5,i)
    semilogy(sort(abs(w_epsilon_increase(:,i)),'descend'))
    xlim([0,1100])
    xlabel('Sorted Index');
    ylabel('Amplitude(dB)');
    title("Excitation Distribution \epsilon = " + epsilon_1(i));
end
%% Finer grid on epsilon over [epsilon_ls_cvx,0.56237], the upper bound is 
% the value whose result start to loose the boundary between zeros and nonzeros

epsilon_2 = linspace(epsilon_ls_cvx,0.56237,25);
w_epsilon_increase2 = zeros(N,size(epsilon_2,2));

%%
for i = 1:1:size(epsilon_2,2)
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S*w - p_d,2) <= epsilon_2(i);
cvx_end
w_epsilon_increase2(:,i) = w;
end

%% plot the l1 norm of the solutions
[N,M] = size(w_epsilon_increase2);
l1_norms2 = zeros(1,M);
for i = 1:M
    l1_norms2(i) = norm(w_epsilon_increase2(:,i),1);
end
figure
plot(epsilon_2(2:end),l1_norms2(2:end));
xlabel('Epsilon');
ylabel('L1 Norm of Excitation Vector');
title('Norm of Vectors in W Against Epsilon');
%%
figure
for i =1:M
    subplot(5,5,i)
    stem(abs(w_epsilon_increase2(:,i)))
    xlim([0,1100])
    xlabel('Elements index');
    ylabel('Amplitude');
    title("Excitation distribution \epsilon = " + epsilon_2(i));
end
%%
figure
for i =1:M
    subplot(5,5,i)
    semilogy(sort(abs(w_epsilon_increase2(:,i)),'descend'))
    xlim([0,1100])
    xlabel('Sorted Index');
    ylabel('Amplitude(dB)');
    title("Excitation Distribution \epsilon = " + epsilon_2(i));
end