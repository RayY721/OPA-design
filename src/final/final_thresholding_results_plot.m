%% Plot for Thresholding (With obtained results, making the plots)
close all
clc
clear
load('LASSO_results_and_epsilons_used.mat')
%% Parameters
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
S = exp(1i*param.k*sin(theta)*x');            % S matr  ix
%% Set the desired beampattern
beam.left = -0.035;
beam.right = 0.035;
p_d = upperboundgen(fov.left,fov.right,[beam.left beam.right],L);
%% Thresholding the result, choosing the one with epsilon equals 0.50558
epsilon_thresholding = epsilon_all(40);
w_2thresholding = w_epsilon_increase_all(:,40);
[W_thresheld, E_thresheld] = threshold(w_2thresholding);

%% Figure out what k value should be used
figure
for i=1:1:11
    subplot(6,2,i)
    plot(angle,20*log10(abs(S*W_thresheld(:,100*(i-1)+1))))
    title("k = " + nnz(W_thresheld(:,100*(i-1)+1)))
end
%% Beampattern at k = 1100, 600, 100
figure
a1 = plot(angle,20*log10(abs(S*W_thresheld(:,1))),'LineWidth', 2);
M1 = "k = " +nnz(W_thresheld(:,1));
hold on
a2 = plot(angle,20*log10(abs(S*W_thresheld(:,601))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "k = " +nnz(W_thresheld(:,601));
a3 = plot(angle,20*log10(abs(S*W_thresheld(:,1001))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M3 = "k = " +nnz(W_thresheld(:,1001));
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampatterns after thresholding to different sparsity level k','fontsize',38);
% grid on;
legend([a1,a2,a3],[M1,M2,M3],'FontSize',39)
xlim([-18 18])
ylim([-70 0])

%% Spatial distribution of excitation at k = 1100, 600, 100
figure
a1 = stem(abs(W_thresheld(:,1)),'LineWidth', 2);
M1 = "k = " +nnz(W_thresheld(:,1));
hold on 
a2 = stem(abs(W_thresheld(:,601)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "k = " +nnz(W_thresheld(:,601));
% a3 = stem(abs(W_thresheld(:,1001)),'Color',[0.4660 0.6740 0.1880],'LineWidth', 3);
% M3 = "k = " +nnz(W_thresheld(:,1001));
set(gca, 'FontSize', 24);
xlabel('Index','fontsize',37);
ylabel('|\mathbf{w}_k|','fontsize',36);
title('The excitation distribution with different values of epsilon','fontsize',36)
legend([a1,a2],[M1,M2],'FontSize',39)
xlim([0,1100])

figure
a1 = stem(abs(W_thresheld(:,1)),'LineWidth', 2);
M1 = "k = " +nnz(W_thresheld(:,1));
hold on 
% a2 = stem(abs(W_thresheld(:,601)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
% M2 = "k = " +nnz(W_thresheld(:,601));
a3 = stem(abs(W_thresheld(:,1001)),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M3 = "k = " +nnz(W_thresheld(:,1001));
set(gca, 'FontSize', 24);
xlabel('Index','fontsize',37);
ylabel('|\mathbf{w}_k|','fontsize',36);
title('The excitation distribution with different values of epsilon','fontsize',36)
legend([a1,a3],[M1,M3],'FontSize',39)
xlim([0,1100])

%% Plot the error distribution k = 100
figure
% a1 = plot(angle,20*log10(abs(S*(W_thresheld(:,1)-W_thresheld(:,601)))),'LineWidth', 2);
a1 = semilogy(angle,abs(S*(W_thresheld(:,1)-W_thresheld(:,1001))),'LineWidth', 2);
hold on 
% M1 = "k = " +nnz(W_thresheld(:,601));
semilogy(angle,mean(abs(S*(W_thresheld(:,1)-W_thresheld(:,1001))))*ones(1,L),'LineWidth',3)
% yline(10*log10(mean(abs(S*(W_thresheld(:,1)-W_thresheld(:,601))))),'Color',[0.9290 0.6940 0.1250],'LineWidth',3)
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('$|\hat{\textbf{\textit{p}}} - \hat{\textbf{\textit{p}}}_{100}|$','Interpreter','latex','fontsize',38);
% ylabel('|\bf\it p ','fontsize',38);
xlim([-18 18])
ylim([1e-4 1])
%% Plot the error distribution k = 500
figure
% a1 = plot(angle,20*log10(abs(S*(W_thresheld(:,1)-W_thresheld(:,601)))),'LineWidth', 2);
a1 = semilogy(angle,abs(S*(W_thresheld(:,1)-W_thresheld(:,601))),'LineWidth', 2);
hold on 
% M1 = "k = " +nnz(W_thresheld(:,601));
semilogy(angle,mean(abs(S*(W_thresheld(:,1)-W_thresheld(:,601))))*ones(1,L),'LineWidth',3)
% yline(10*log10(mean(abs(S*(W_thresheld(:,1)-W_thresheld(:,601))))),'Color',[0.9290 0.6940 0.1250],'LineWidth',3)
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('$|\hat{\textbf{\textit{p}}} - \hat{\textbf{\textit{p}}}_{500}|$','Interpreter','latex','fontsize',38);
% ylabel('|\bf\it p ','fontsize',38);
xlim([-18 18])
ylim([1e-4 1])
% %%
% figure
% a2 = plot(angle,20*log10(abs(S*(W_thresheld(:,1)-W_thresheld(:,1001)))),'LineWidth', 2);
% % M2 = "k = " +nnz(W_thresheld(:,1001));
% yline(20*log10(mean(abs(S*(W_thresheld(:,1)-W_thresheld(:,1001))))),'Color',[0.9290 0.6940 0.1250],'LineWidth',3)
% set(gca, 'FontSize', 26);
% xlabel('\theta','fontsize',39);
% ylabel('Intensity(dB)','fontsize',38);
% % title('Beampatterns after thresholding to different sparsity level k','fontsize',38);
% % grid on;
% % legend(a1,M1,'FontSize',39)
% xlim([-18 18])
% ylim([-70 0])

%% Plot the entries have been set to 0
figure
stem(abs(W_thresheld(:,1)-W_thresheld(:,601)),'LineWidth',2);
set(gca, 'FontSize', 26);
xlabel('Index','fontsize',39);
ylabel('$|\hat{\textbf{\textit{w}}} - \hat{\textbf{\textit{w}}}_{500}|$','Interpreter','latex','fontsize',38);
xlim([0,1100])
ylim([0,2.5e-3])
figure
stem(abs(W_thresheld(:,1)-W_thresheld(:,1001)),'LineWidth',2);
set(gca, 'FontSize', 26);
xlabel('Index','fontsize',39);
ylabel('$|\hat{\textbf{\textit{w}}} - \hat{\textbf{\textit{w}}}_{100}|$','Interpreter','latex','fontsize',38);
xlim([0,1100])
ylim([0,2.5e-3])

%% Mesh plot for S^HS
% Display the S^H * S 
[X,Y] = meshgrid(1:N,1:N);
figure
surf(X,Y,real(S'*S))
grid off
set(gca, 'FontSize', 26);
xlabel('Index of Row','FontSize',39)
ylabel('Index of Column','FontSize',39)
zlabel('Amplitude','FontSize',39)

%% Plot the beampattern with expected k = 500
figure
a2 = plot(angle,20*log10(abs(S*W_thresheld(:,601))),'LineWidth', 2);
M2 = "k = " +nnz(W_thresheld(:,601));
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampatterns after thresholding to different sparsity level k','fontsize',38);
yline(20*log10(mean(abs(S*(W_thresheld(:,1)-W_thresheld(:,601))))),'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
yline(20*log10(norm(W_thresheld(:,1)-W_thresheld(:,601))),'LineWidth',2)

%% Plot the beampattern with expected k = 10000
figure
a2 = plot(angle,20*log10(abs(S*W_thresheld(:,1001))),'LineWidth', 2);
M2 = "k = " +nnz(W_thresheld(:,1001));
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampatterns after thresholding to different sparsity level k','fontsize',38);
yline(20*log10(mean(abs(S*(W_thresheld(:,1)-W_thresheld(:,1001))))),'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
yline(20*log10(norm(W_thresheld(:,1)-W_thresheld(:,1001))),'LineWidth',2)