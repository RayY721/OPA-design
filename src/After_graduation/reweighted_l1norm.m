%% This script implements the reweighted l1 norm minimization 
% The reweigthed l1 norm minimization outperforms the ordinary l1 norm
% minimization. By controlling the hyperparameter, the solution shows
% different level of sparsity
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
%% Stored results and parameters
% Instead of running the algorithm, results evaluated with according
% parameters are stored in the following .mat file
% loading this .mat file can skip running the algo and jumping into the
% section where the results are examed. 
load('reweighted_results.mat')
%% Reweighted l1 norm minimization
% Use a order 3 tensor to store the results of reweighted algo, 
% For single epsilon value, the result is a NxP matrix 
% For J = size(epsilon_reweighted,2) number of epsilon values, the result
% should be an order 3 tensor: (N,P,J) 

% An vector (epsilon_reweighted) need to be created to specify the values
% of epsilon. 
% An example is as following: (the sections for exam the results assume
% this vector has four entries) 
epsilon_reweighted = [0.394671089074041	0.505582211203687	1.03088596408733	1.37162160700801];
J = size(epsilon_reweighted,2);             % Number of epsilon values
P = 5;                                      % Algo stops after 5 iterations
W_reweighted = zeros(N,P,J);                % An order 3 tensor to store the result
optimal_value = zeros(P,J);

for j =1:1:J
    Z = eye(N);
    for i = 1:1:P                           
    cvx_begin
        variable w(N) complex
        minimize(norm(Z*w,1))               % The norm is weighted by matrix Z
        subject to
            norm(S*w - p_d,2) <= epsilon_reweighted(j);       
    cvx_end
    W_reweighted(:,i,j) = w;
    Z = inv(diag(abs(w) + 0.0002));
    optimal_value(i,j) = cvx_optval;        % The optimal value of this iteration
    end
end

%% Exam the result with epsilon = 0.5 (j = 2)

% Plot the resulting beampattern at epsilon = 0.5. This plot demonstrate
% the effect of iteratively weighting the l1 norm
figure
a1 = plot(angle,20*log10(abs(S*W_reweighted(:,1,2))),'LineWidth', 2);
M1 = "No.1 iteration";
hold on
a2 = plot(angle,20*log10(abs(S*W_reweighted(:,5,2))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "No.5 iteration";
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampatterns at the First and the Last Iteration','fontsize',38);
legend([a1,a2],[M1,M2],'FontSize',39)
xlim([-18 18])
ylim([-70 0])

%% Exam the results with different epsilon values
% Plot the sorted amplitude distribution for different epsilon
figure
a1 = semilogy(sort(abs(W_reweighted(:,5,1)),'descend'),'LineWidth', 2);
M1 = "\epsilon = " + epsilon_reweighted(1);
hold on
a2 = semilogy(sort(abs(W_reweighted(:,5,2)),'descend'),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "\epsilon = " + epsilon_reweighted(2);
a3 = semilogy(sort(abs(W_reweighted(:,5,3)),'descend'),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M3 = "\epsilon = " + epsilon_reweighted(3);
a4 = semilogy(sort(abs(W_reweighted(:,5,4)),'descend'),'LineWidth', 2);
M4 = "\epsilon = " + epsilon_reweighted(4);
set(gca, 'FontSize', 26);
xlabel('Index','fontsize',39);
ylabel('Magnitude(dB)','fontsize',38);
title("The Sorted Amplitude Distribution",'fontsize',38)
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([0,1100])

% plot the spatial excitation distribution in seperate plot 
figure
subplot(2,1,1)
a1 = stem(abs(W_reweighted(:,5,1)),'LineWidth', 1);
M1 = "\epsilon = " + epsilon_reweighted(1);
hold on 
a2 = stem(abs(W_reweighted(:,5,2)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 1);
M2 = "\epsilon = " + epsilon_reweighted(2);
set(gca, 'FontSize', 24);
xlabel('Index','fontsize',37);
ylabel('Amplitude','fontsize',36);
% title('The excitation distribution with different values of epsilon','fontsize',36)
legend([a1,a2],[M1,M2],'FontSize',39)
xlim([0,1100])
ylim([0,4e-3])
subplot(2,1,2)
a3 = stem(abs(W_reweighted(:,5,3)),'Color',[0.4660 0.6740 0.1880],'LineWidth', 1);
M3 = "\epsilon = " + epsilon_reweighted(3);
hold on
a4 = stem(abs(W_reweighted(:,5,4)),'LineWidth', 1);
M4 = "\epsilon = " + epsilon_reweighted(4);
set(gca, 'FontSize', 24);
xlabel('Index','fontsize',37);
ylabel('Amplitude','fontsize',36);
% title('The excitation distribution with different values of epsilon','fontsize',36)
legend([a3,a4],[M3,M4],'FontSize',39)
xlim([0,1100])
ylim([0,4e-3])

% Plot all beampatterns in one figure
figure;
a1 = plot(angle,20*log10(abs(S*W_reweighted(:,5,1))),'LineWidth', 2);
M1 = "\epsilon = " + epsilon_reweighted(1);
hold on 
a2 = plot(angle,20*log10(abs(S*W_reweighted(:,5,2))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "\epsilon = " + epsilon_reweighted(2);
a3 = plot(angle,20*log10(abs(S*W_reweighted(:,5,3))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M3 = "\epsilon = " + epsilon_reweighted(3);
a4 = plot(angle,20*log10(abs(S*W_reweighted(:,5,4))),'LineWidth', 3);
M4 = "\epsilon = " + epsilon_reweighted(4);
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampattern','fontsize',38);
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([-18 18])
ylim([-70 0])

% The investigation of how epsilon is used on the beampattern matching
figure
a1 = stem(angle,100*abs(S*W_reweighted(:,5,1)-p_d).^2./epsilon_reweighted(1)^2,'LineWidth', 3);
M1 = "\epsilon = " + epsilon_reweighted(1);
hold on
a2 = stem(angle,100*abs(S*W_reweighted(:,5,2)-p_d).^2./epsilon_reweighted(2)^2,'Color',[0.9290 0.6940 0.1250],'LineWidth', 3);
M2 = "\epsilon = " + epsilon_reweighted(2);
a3 = stem(angle,100*abs(S*W_reweighted(:,5,3)-p_d).^2./epsilon_reweighted(3)^2,'Color',[0.4660 0.6740 0.1880],'LineWidth', 3);
M3 = "\epsilon = " + epsilon_reweighted(3);
a4 = stem(angle,100*abs(S*W_reweighted(:,5,4)-p_d).^2./epsilon_reweighted(4)^2,'LineWidth', 3);
M4 = "\epsilon = " + epsilon_reweighted(4);
set(gca, 'FontSize', 24);
xlabel('\theta','fontsize',37);
ylabel('The percentage (%)','fontsize',36);
% title('The normalized square error (\epsilon^2) distribution with different \epsilon','fontsize',36)
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([-1,1])
ylim([0,50])
