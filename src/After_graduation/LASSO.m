%% This script evaluate the epsilon at different values
close all
clc
clear
%% Parameters
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;

N = 1100;               % to have a sufficient aperture coverage
d = param.lambda;       % A grid with one wavelength spacing
L = 1100;               % so that the sampling interval is small enough
% set the range on which the beampattern matches
fov.left = -18;
fov.right = 18;

% generating the points on which (the grid used for matching) 
angle = linspace(fov.left,fov.right,L)';  % in degree
theta = angle*pi/180;

% The grid to discretize the aperture
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix

% Generating the desired beampattern
% The beamwidth is set to be 0.07 degree
beam.left = -0.035;
beam.right = 0.035;
% Call a function to generate a desred beampattern vector 
p_d = desiredbeam(fov.left,fov.right,[beam.left beam.right],L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least square solution investigation
%% LS solution with pinv()

% Calculate the LS solution with pinv()
w_ls = pinv(S)*p_d;     % This solution gives the LS error
epsilon_ls = norm(S*w_ls - p_d,2);
epsilon_max = norm(p_d,2);

% plot the beampattern
figure
plot(angle,20*log10(abs(S*w_ls)));
set(gca,'FontSize',12);
xlabel('\theta','FontSize',16);
ylabel('Intensity(dB)','FontSize',16);
title('Beampattern','FontSize',16);
% plot the distribution
figure
stem(abs(w_ls))
xlim([0,1100])
set(gca,'FontSize',12);
xlabel('Index of grid points','FontSize',16);
ylabel('Amplitude','FontSize',16);
title('Excitation over the Grid','FontSize',16);

%% LS solution with CVX

% CVX for LS solution
cvx_begin
    variable w_ls_cvx(N) complex
    minimize(norm(S*w_ls_cvx - p_d))
cvx_end
epsilon_ls_cvx = norm(S*w_ls_cvx - p_d,2);

% plot the beampattern
figure
plot(angle,20*log10(abs(S*w_ls_cvx)));
set(gca,'FontSize',12);
xlabel('\theta','FontSize',16);
ylabel('Intensity(dB)','FontSize',16);
title('Beampattern','FontSize',16);

% plot the distribution
figure
stem(abs(w_ls_cvx))
xlim([0,1100])
set(gca,'FontSize',12);
xlabel('Index of grid points','FontSize',16);
ylabel('Amplitude','FontSize',16);
title('Excitation over the Grid','FontSize',16);


%% LASSO problem
load('LASSO_parameters.mat')    % load the epsilon values for later use
% It also includes the solution (w) for different epsilon values, so no
% need to run the following section unless you want to change something

%%
w_epsilon_increase = zeros(N,size(epsilon_all,2));
for i = 1:1:size(epsilon_all,2)
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S*w - p_d,2) <= epsilon_all(i);
cvx_end
w_epsilon_increase(:,i) = w;
end

%% Exam the results

% plot the epsilon values
figure
plot(epsilon_all,'o')
grid on
set(gca,'FontSize',12);
xlabel('Index','FontSize',16);
ylabel('\epsilon','FontSize',16);
% title('The Plot of the Series \epsilon Values','FontSize',16);
yline(max(epsilon_all),'LineWidth',1,'Color','r')
yline(min(epsilon_all),'LineWidth',1,'Color','r')

% plot the l1 norm of the solutions
[N,M] = size(w_epsilon_increase);
l1_norms = zeros(1,M);
for i = 1:M
    l1_norms(i) = norm(w_epsilon_increase(:,i),1);
end
figure
semilogy(epsilon_all,l1_norms,'LineWidth',2);
set(gca,'FontSize',12);
xlabel('\epsilon','FontSize',16);
ylabel('$||\textbf{w}||_1$','Interpreter','latex','FontSize',16);
% title('Norm of Vectors in W Against Epsilon','FontSize',16);
xlim([min(epsilon_all),max(epsilon_all)])

% select a few epsilon value to have a closer look
selection_vector = [9 40 59 67];

% Plot all beampattern in one figure
figure;
a1 = plot(angle,20*log10(abs(S*w_epsilon_increase(:,selection_vector(1)))),'LineWidth', 2);
M1 = "\epsilon = "+epsilon_all(selection_vector(1));
hold on 
a2 = plot(angle,20*log10(abs(S*w_epsilon_increase(:,selection_vector(2)))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "\epsilon = "+epsilon_all(selection_vector(2));
a3 = plot(angle,20*log10(abs(S*w_epsilon_increase(:,selection_vector(3)))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M3 = "\epsilon = "+epsilon_all(selection_vector(3));
a4 = plot(angle,20*log10(abs(S*w_epsilon_increase(:,selection_vector(4)))),'LineWidth', 3);
M4 = "\epsilon = "+epsilon_all(selection_vector(4));
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
% title('Beampattern of results with different values of epsilon','fontsize',38);
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([-18 18])
ylim([-70 0])

% Plot all excitation distribution in one figure
figure
a1 = stem(abs(w_epsilon_increase(:,selection_vector(1))),'LineWidth', 1);
M1 = "\epsilon = "+epsilon_all(selection_vector(1));
hold on 
a2 = stem(abs(w_epsilon_increase(:,selection_vector(2))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 1);
M2 = "\epsilon = "+epsilon_all(selection_vector(2));
a3 = stem(abs(w_epsilon_increase(:,selection_vector(3))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 1);
M3 = "\epsilon = "+epsilon_all(selection_vector(3));
a4 = stem(abs(w_epsilon_increase(:,selection_vector(4))),'LineWidth', 1);
M4 = "\epsilon = "+epsilon_all(selection_vector(4));
set(gca, 'FontSize', 24);
xlabel('Index','fontsize',37);
ylabel('Amplitude','fontsize',36);
% title('The excitation distribution with different values of epsilon','fontsize',36)
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([0,1100])
ylim([0,4e-3])


% Plot all sorted amplitude distribution in one figure
figure
a1 = semilogy(sort(abs(w_epsilon_increase(:,selection_vector(1))),'descend'),'LineWidth', 2);
M1 = "\epsilon = "+epsilon_all(selection_vector(1));
hold on
a2 = semilogy(sort(abs(w_epsilon_increase(:,selection_vector(2))),'descend'),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "\epsilon = "+epsilon_all(selection_vector(2));
a3 = semilogy(sort(abs(w_epsilon_increase(:,selection_vector(3))),'descend'),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M3 = "\epsilon = "+epsilon_all(selection_vector(3));
a4 = semilogy(sort(abs(w_epsilon_increase(:,selection_vector(4))),'descend'),'LineWidth', 2);
M4 = "\epsilon = "+epsilon_all(selection_vector(4));
set(gca, 'FontSize', 24);
xlabel('Index','fontsize',37);
ylabel('Magnitude(dB)','fontsize',36);
title("The Sorted Amplitude Distribution",'fontsize',38)
% title('The sorted unit power excitation distribution with different values of epsilon','fontsize',36)
grid on
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([0,1100])

% Plot all normalized squared epsilon distribution is one figure
figure
a1 = stem(angle,100*abs(S*w_epsilon_increase(:,selection_vector(1))-p_d).^2./epsilon_all(selection_vector(1))^2,'LineWidth', 3);
M1 = "\epsilon = "+epsilon_all(selection_vector(1));
hold on
a2 = stem(angle,100*abs(S*w_epsilon_increase(:,selection_vector(2))-p_d).^2./epsilon_all(selection_vector(2))^2,'Color',[0.9290 0.6940 0.1250],'LineWidth', 3);
M2 = "\epsilon = "+epsilon_all(selection_vector(2));
a3 = stem(angle,100*abs(S*w_epsilon_increase(:,selection_vector(3))-p_d).^2./epsilon_all(selection_vector(3))^2,'Color',[0.4660 0.6740 0.1880],'LineWidth', 3);
M3 = "\epsilon = "+epsilon_all(selection_vector(3));
a4 = stem(angle,100*abs(S*w_epsilon_increase(:,selection_vector(4))-p_d).^2./epsilon_all(selection_vector(4))^2,'LineWidth', 3);
M4 = "\epsilon = "+epsilon_all(selection_vector(4));
set(gca, 'FontSize', 24);
xlabel('\theta','fontsize',37);
ylabel('The percentage (%)','fontsize',36);
% title('The normalized square error (\epsilon^2) distribution with different \epsilon','fontsize',36)
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([-1,1])


%% Plot the least square solution for comparison

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
