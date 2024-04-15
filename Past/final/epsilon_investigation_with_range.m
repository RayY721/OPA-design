%% This script exam the choice of the epsilon
close all
clc
clear
load('epsilon_investigation_with_range.mat')
epsilon_ = [0.4 0.65 0.9 1.15 1.4];

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
S = exp(1i*param.k*sin(theta)*x');            % S matrix

%% Set the desired beampattern
beam.left = -0.035;
beam.right = 0.035;
p_d = upperboundgen(fov.left,fov.right,[beam.left beam.right],L);
p_d = p_d*1i;
% %% Modify the S matrix and reference pattern
% p_d_modified = p_d;
% S = S;
%% Determine the range of the epsilon
w_ls = pinv(S)*p_d;     % This solution gives the LS error
epsilon_ls = norm(S*w_ls - p_d,2);
epsilon_max = norm(p_d,2);
%% CVX for LS solution
cvx_begin
    variable w_ls_cvx(N) complex
    minimize(norm(S*w_ls_cvx - p_d))
cvx_end
%% 
% epsilon_ = linspace(epsilon_ls,epsilon_max,5);
w_epsilon_increase = zeros(N,size(epsilon_,2));
for i = 1:1:size(epsilon_,2)
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S*w - p_d,2) <= epsilon_(i);
cvx_end
w_epsilon_increase(:,i) = w;
end
% w_epsilon_increase(:,1) = w_ls; 
%%
figure
for i = 1:1:size(epsilon_,2)
    subplot(2,3,i)
    plot(angle,20*log10(abs(S*w_epsilon_increase(:,i))))
    xlabel('Angle in degree')
    ylabel('Intensity in dB')
    % ylim([-60,0])
    xlim([-18,18])
    title("relaxation of epsilon =" + epsilon_(i))
end
%% plot four beampatterns in one 
figure;
% plot(angle,20*log10(abs(S*w)))
a1 = plot(angle,20*log10(abs(S*w_epsilon_increase(:,1))),'LineWidth', 2);
M1 = "\epsilon = "+epsilon_(1);
hold on 
a2 = plot(angle,20*log10(abs(S*w_epsilon_increase(:,2))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "\epsilon = "+epsilon_(2);
a3 = plot(angle,20*log10(abs(S*w_epsilon_increase(:,3))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M3 = "\epsilon = "+epsilon_(3);
a4 = plot(angle,20*log10(abs(S*w_epsilon_increase(:,4))),'LineWidth', 2);
M4 = "\epsilon = "+epsilon_(4);
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampattern of results with different values of epsilon','fontsize',38);
% grid on;
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([-18 18])
ylim([-60 0])
%% plot the normalized epsilon distribution
figure
a1 = stem(angle,100*abs(S*w_epsilon_increase(:,1)-p_d).^2./epsilon_(1)^2,'LineWidth', 2);
M1 = "\epsilon = "+epsilon_(1);
hold on
a2 = stem(angle,100*abs(S*w_epsilon_increase(:,2)-p_d).^2./epsilon_(2)^2,'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "\epsilon = "+epsilon_(2);
a3 = stem(angle,100*abs(S*w_epsilon_increase(:,3)-p_d).^2./epsilon_(3)^2,'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M3 = "\epsilon = "+epsilon_(3);
a4 = stem(angle,100*abs(S*w_epsilon_increase(:,4)-p_d).^2./epsilon_(4)^2,'LineWidth', 2);
M4 = "\epsilon = "+epsilon_(4);
set(gca, 'FontSize', 24);
xlabel('\theta','fontsize',37);
ylabel('The percentage of \epsilon^2 (%)','fontsize',36);
title('The normalized square error (\epsilon^2) distribution with different \epsilon','fontsize',36)
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([-18,18])
figure
a1 = stem(angle,100*abs(S*w_epsilon_increase(:,1)-p_d).^2./epsilon_(1)^2,'LineWidth', 3);
M1 = "\epsilon = "+epsilon_(1);
hold on
a2 = stem(angle,100*abs(S*w_epsilon_increase(:,2)-p_d).^2./epsilon_(2)^2,'Color',[0.9290 0.6940 0.1250],'LineWidth', 3);
M2 = "\epsilon = "+epsilon_(2);
a3 = stem(angle,100*abs(S*w_epsilon_increase(:,3)-p_d).^2./epsilon_(3)^2,'Color',[0.4660 0.6740 0.1880],'LineWidth', 3);
M3 = "\epsilon = "+epsilon_(3);
a4 = stem(angle,100*abs(S*w_epsilon_increase(:,4)-p_d).^2./epsilon_(4)^2,'LineWidth', 3);
M4 = "\epsilon = "+epsilon_(4);
set(gca, 'FontSize', 24);
xlabel('\theta','fontsize',37);
ylabel('The percentage of \epsilon^2 (%)','fontsize',36);
title('The normalized square error (\epsilon^2) distribution with different \epsilon','fontsize',36)
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([-0.8,0.8])
%% The exciation part
% figure
% a1 = stem(abs(w_epsilon_increase(:,1)),'LineWidth', 1);
% M1 = "\epsilon = "+epsilon_(1);
% hold on
% a2 = stem(abs(w_epsilon_increase(:,2)),'LineWidth', 1);
% M2 = "\epsilon = "+epsilon_(2);
% a3 = stem(abs(w_epsilon_increase(:,3)),'LineWidth', 1);
% M3 = "\epsilon = "+epsilon_(3);
% a4 = stem(abs(w_epsilon_increase(:,4)),'LineWidth', 1);
% M4 = "\epsilon = "+epsilon_(4);
% xlabel('\theta','fontsize',21);
% ylabel('The square error at sampling points','fontsize',20);
% title('The normalized error distribution with different value of epsilon','fontsize',20)
% legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',21)
% xlim([0,1100])

figure
a1 = stem(abs(w_epsilon_increase(:,1)),'LineWidth', 1);
M1 = "\epsilon = "+epsilon_(1);
hold on 
a2 = stem(abs(w_epsilon_increase(:,2)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 1);
M2 = "\epsilon = "+epsilon_(2);
a3 = stem(abs(w_epsilon_increase(:,3)),'Color',[0.4660 0.6740 0.1880],'LineWidth', 1);
M3 = "\epsilon = "+epsilon_(3);
a4 = stem(abs(w_epsilon_increase(:,4)),'LineWidth', 1);
M4 = "\epsilon = "+epsilon_(4);
set(gca, 'FontSize', 24);
xlabel('The index of elements','fontsize',37);
ylabel('The amplitude of excitation','fontsize',36);
title('The excitation distribution with different values of epsilon','fontsize',36)
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([0,1100])

%% Unit power excitation 
figure
a1 = stem(abs(w_epsilon_increase(:,1)./norm(w_epsilon_increase(:,1),2)),'LineWidth', 1);
M1 = "\epsilon = "+epsilon_(1);
hold on 
a2 = stem(abs(w_epsilon_increase(:,2))./norm(w_epsilon_increase(:,2),2),'Color',[0.9290 0.6940 0.1250],'LineWidth', 1);
M2 = "\epsilon = "+epsilon_(2);
a3 = stem(abs(w_epsilon_increase(:,3))./norm(w_epsilon_increase(:,3),2),'Color',[0.4660 0.6740 0.1880],'LineWidth', 1);
M3 = "\epsilon = "+epsilon_(3);
a4 = stem(abs(w_epsilon_increase(:,4))./norm(w_epsilon_increase(:,4),2),'LineWidth', 1);
M4 = "\epsilon = "+epsilon_(4);
set(gca, 'FontSize', 24);
xlabel('The index of elements','fontsize',37);
ylabel('The amplitude of excitation','fontsize',36);
title('The unit power excitation distribution with different values of epsilon','fontsize',36)
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([0,1100])

%% The sorted distribution of excitation
figure
a1 = semilogy(sort(abs(w_epsilon_increase(:,1)),'descend')./norm(w_epsilon_increase(:,1),2),'LineWidth', 1);
% a1 = semilogy(sort(abs(w_epsilon_increase(:,1)),'descend'),'LineWidth', 1);
M1 = "\epsilon = "+epsilon_(1);
hold on
a2 = semilogy(sort(abs(w_epsilon_increase(:,2)),'descend')./norm(w_epsilon_increase(:,2),2),'Color',[0.9290 0.6940 0.1250],'LineWidth', 1);
% a2 = semilogy(sort(abs(w_epsilon_increase(:,2)),'descend'),'Color',[0.9290 0.6940 0.1250],'LineWidth', 1);
M2 = "\epsilon = "+epsilon_(2);
a3 = semilogy(sort(abs(w_epsilon_increase(:,3)),'descend')./norm(w_epsilon_increase(:,3),2),'Color',[0.4660 0.6740 0.1880],'LineWidth', 1);
% a3 = semilogy(sort(abs(w_epsilon_increase(:,3)),'descend'),'Color',[0.4660 0.6740 0.1880],'LineWidth', 1);
M3 = "\epsilon = "+epsilon_(3);
a4 = semilogy(sort(abs(w_epsilon_increase(:,4)),'descend')./norm(w_epsilon_increase(:,4),2),'LineWidth', 1);
% a4 = semilogy(sort(abs(w_epsilon_increase(:,4)),'descend'),'LineWidth', 1);
M4 = "\epsilon = "+epsilon_(4);
set(gca, 'FontSize', 24);
xlabel('The index of elements','fontsize',37);
ylabel('The amplitude of excitation(dB)','fontsize',36);
title('The sorted unit power excitation distribution with different values of epsilon','fontsize',36)
grid on
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([0,1100])

%% Beampattern with unit power
figure
for i = 1:1:size(epsilon_,2)
    subplot(2,3,i)
    plot(angle,20*log10(abs(S*w_epsilon_increase(:,i)./norm(w_epsilon_increase(:,i),2))))
    xlabel('Angle in degree')
    ylabel('Intensity in dB')
    % ylim([-60,0])
    xlim([-18,18])
    title("relaxation of epsilon =" + epsilon_(i))
end