%% With obtained results, making the plots
%% This script exam the choice of the epsilon
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
S = exp(1i*param.k*sin(theta)*x');            % S matrix
%% Set the desired beampattern
beam.left = -0.035;
beam.right = 0.035;
p_d = upperboundgen(fov.left,fov.right,[beam.left beam.right],L);
%% plot the seires of epsilon
figure
plot(epsilon_all,'o')
grid on
set(gca,'FontSize',12);
xlabel('Index','FontSize',16);
ylabel('\epsilon','FontSize',16);
% title('The Plot of the Series \epsilon Values','FontSize',16);
yline(max(epsilon_all),'LineWidth',1,'Color','r')
yline(min(epsilon_all),'LineWidth',1,'Color','r')

%% plot the l1 norm of the solutions
[N,M] = size(w_epsilon_increase_all);
l1_norms = zeros(1,M);
for i = 1:M
    l1_norms(i) = norm(w_epsilon_increase_all(:,i),1);
end
figure
semilogy(epsilon_all,l1_norms,'LineWidth',2);
set(gca,'FontSize',12);
xlabel('\epsilon','FontSize',16);
ylabel('$||\textbf{w}||_1$','Interpreter','latex','FontSize',16);
% title('Norm of Vectors in W Against Epsilon','FontSize',16);
xlim([min(epsilon_all),max(epsilon_all)])

% %% plot all solutions
% figure
% for i =1:M
%     subplot(9,8,i)
%     stem(abs(w_epsilon_increase_all(:,i)))
%     xlim([0,1100])
%     xlabel('Elements index');
%     ylabel('Amplitude');
%     title("Excitation distribution \epsilon = " + epsilon_all(i));
% end
% %%
% figure
% for i =1:M
%     subplot(9,8,i)
%     semilogy(sort(abs(w_epsilon_increase_all(:,i)),'descend'))
%     xlim([0,1100])
%     xlabel('Sorted Index');
%     ylabel('Amplitude(dB)');
%     title("Excitation Distribution \epsilon = " + epsilon_all(i));
% end

%% selected epsilon values
selection_vector = [9 40 59 67];
figure
for i = 1:size(selection_vector,2)
    subplot(2,2,i)
    stem(abs(w_epsilon_increase_all(:,selection_vector(i))))
    xlim([0,1100])
    xlabel('Elements index');
    ylabel('Amplitude');
    title("Excitation distribution \epsilon = " + epsilon_all(selection_vector(i)));
end

figure
for i = 1:size(selection_vector,2)
    subplot(2,2,i)
    semilogy(sort(abs(w_epsilon_increase_all(:,selection_vector(i))),'descend'))
    xlim([0,1100])
    xlabel('Elements index');
    ylabel('Amplitude');
    title("Excitation distribution \epsilon = " + epsilon_all(selection_vector(i)));
end

%% The beampattern in one 
figure;
% plot(angle,20*log10(abs(S*w)))
a1 = plot(angle,20*log10(abs(S*w_epsilon_increase_all(:,selection_vector(1)))),'LineWidth', 2);
M1 = "\epsilon = "+epsilon_all(selection_vector(1));
hold on 
a2 = plot(angle,20*log10(abs(S*w_epsilon_increase_all(:,selection_vector(2)))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "\epsilon = "+epsilon_all(selection_vector(2));
a3 = plot(angle,20*log10(abs(S*w_epsilon_increase_all(:,selection_vector(3)))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M3 = "\epsilon = "+epsilon_all(selection_vector(3));
a4 = plot(angle,20*log10(abs(S*w_epsilon_increase_all(:,selection_vector(4)))),'LineWidth', 3);
M4 = "\epsilon = "+epsilon_all(selection_vector(4));
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
% title('Beampattern of results with different values of epsilon','fontsize',38);
% grid on;
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([-18 18])
ylim([-70 0])

%% The excitation distribution in one
figure
a1 = stem(abs(w_epsilon_increase_all(:,selection_vector(1))),'LineWidth', 1);
M1 = "\epsilon = "+epsilon_all(selection_vector(1));
hold on 
a2 = stem(abs(w_epsilon_increase_all(:,selection_vector(2))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 1);
M2 = "\epsilon = "+epsilon_all(selection_vector(2));
a3 = stem(abs(w_epsilon_increase_all(:,selection_vector(3))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 1);
M3 = "\epsilon = "+epsilon_all(selection_vector(3));
a4 = stem(abs(w_epsilon_increase_all(:,selection_vector(4))),'LineWidth', 1);
M4 = "\epsilon = "+epsilon_all(selection_vector(4));
set(gca, 'FontSize', 24);
xlabel('Index','fontsize',37);
ylabel('Amplitude','fontsize',36);
% title('The excitation distribution with different values of epsilon','fontsize',36)
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([0,1100])
ylim([0,4e-3])

%% The sorted amplitude distribution in one
figure
a1 = semilogy(sort(abs(w_epsilon_increase_all(:,selection_vector(1))),'descend'),'LineWidth', 2);
% a1 = semilogy(sort(abs(w_epsilon_increase(:,1)),'descend'),'LineWidth', 1);
M1 = "\epsilon = "+epsilon_all(selection_vector(1));
hold on
a2 = semilogy(sort(abs(w_epsilon_increase_all(:,selection_vector(2))),'descend'),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
% a2 = semilogy(sort(abs(w_epsilon_increase(:,2)),'descend'),'Color',[0.9290 0.6940 0.1250],'LineWidth', 1);
M2 = "\epsilon = "+epsilon_all(selection_vector(2));
a3 = semilogy(sort(abs(w_epsilon_increase_all(:,selection_vector(3))),'descend'),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
% a3 = semilogy(sort(abs(w_epsilon_increase(:,3)),'descend'),'Color',[0.4660 0.6740 0.1880],'LineWidth', 1);
M3 = "\epsilon = "+epsilon_all(selection_vector(3));
a4 = semilogy(sort(abs(w_epsilon_increase_all(:,selection_vector(4))),'descend'),'LineWidth', 2);
% a4 = semilogy(sort(abs(w_epsilon_increase(:,4)),'descend'),'LineWidth', 1);
M4 = "\epsilon = "+epsilon_all(selection_vector(4));
set(gca, 'FontSize', 24);
xlabel('Index','fontsize',37);
ylabel('Magnitude(dB)','fontsize',36);
title("The Sorted Amplitude Distribution",'fontsize',38)
% title('The sorted unit power excitation distribution with different values of epsilon','fontsize',36)
grid on
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([0,1100])

%% plot the normalized squared epsilon distribution
figure
a1 = stem(angle,100*abs(S*w_epsilon_increase_all(:,selection_vector(1))-p_d).^2./epsilon_all(selection_vector(1))^2,'LineWidth', 3);
M1 = "\epsilon = "+epsilon_all(selection_vector(1));
hold on
a2 = stem(angle,100*abs(S*w_epsilon_increase_all(:,selection_vector(2))-p_d).^2./epsilon_all(selection_vector(2))^2,'Color',[0.9290 0.6940 0.1250],'LineWidth', 3);
M2 = "\epsilon = "+epsilon_all(selection_vector(2));
a3 = stem(angle,100*abs(S*w_epsilon_increase_all(:,selection_vector(3))-p_d).^2./epsilon_all(selection_vector(3))^2,'Color',[0.4660 0.6740 0.1880],'LineWidth', 3);
M3 = "\epsilon = "+epsilon_all(selection_vector(3));
a4 = stem(angle,100*abs(S*w_epsilon_increase_all(:,selection_vector(4))-p_d).^2./epsilon_all(selection_vector(4))^2,'LineWidth', 3);
M4 = "\epsilon = "+epsilon_all(selection_vector(4));
set(gca, 'FontSize', 24);
xlabel('\theta','fontsize',37);
ylabel('The percentage (%)','fontsize',36);
% title('The normalized square error (\epsilon^2) distribution with different \epsilon','fontsize',36)
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([-1,1])

%% plot the LS solution

figure
stem(abs(w_ls_cvx),'LineWidth', 2)
set(gca,'FontSize', 26);
xlim([0,1100])
xlabel('Index','fontsize',39);
ylabel('The Excitation Amplitude','fontsize',38);
title('The Spatial Amplitude Distribution with \epsilon_{LS}','fontsize',38);

figure
semilogy(sort(abs(w_ls_cvx),'descend'),'LineWidth', 2)
set(gca,'FontSize', 26);
xlim([0,1100])
xlabel('Index','fontsize',39);
ylabel('The Excitation Amplitude(dB)','fontsize',38);
title('The Sorted Amplitude Distribution with \epsilon_{LS}','fontsize',38);