%% Get the readings
close all
clc
clear

% Get the readings, use a larger L
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;
N = 1100;                % to have a aperture of 1020 lambda
d = param.lambda;
L = 200001;             % The resolution is at least 0.035 degree
fov.left = -18;
fov.right = 18;
Res = (fov.right - fov.left)/(L - 1);
angle = linspace(fov.left,fov.right,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix


%% Reweighted results
load('reweighted_results.mat')

figure;
a1 = plot(angle,20*log10(abs(S*W_reweighted(:,5,1))),'LineWidth', 2);
M1 = "\epsilon = " + epsilon_reweighted(1);
% hold on 
figure
a2 = plot(angle,20*log10(abs(S*W_reweighted(:,5,2))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M1 = "\epsilon = " + epsilon_reweighted(2);
figure
a3 = plot(angle,20*log10(abs(S*W_reweighted(:,5,3))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M1 = "\epsilon = " + epsilon_reweighted(3);
figure
a4 = plot(angle,20*log10(abs(S*W_reweighted(:,5,4))),'LineWidth', 3);
M1 = "\epsilon = " + epsilon_reweighted(4);
%%
disp("l1 norm of result when epsilon = "+ epsilon_reweighted(1)+ " is " + norm(W_reweighted(:,5,1),1))
disp("l2 norm of result when epsilon = "+ epsilon_reweighted(1)+ " is " + norm(W_reweighted(:,5,1),2))

disp("l1 norm of result when epsilon = "+ epsilon_reweighted(2)+ " is " + norm(W_reweighted(:,5,2),1))
disp("l2 norm of result when epsilon = "+ epsilon_reweighted(2)+ " is " + norm(W_reweighted(:,5,2),2))

disp("l1 norm of result when epsilon = "+ epsilon_reweighted(3)+ " is " + norm(W_reweighted(:,5,3),1))
disp("l2 norm of result when epsilon = "+ epsilon_reweighted(3)+ " is " + norm(W_reweighted(:,5,3),2))

disp("l1 norm of result when epsilon = "+ epsilon_reweighted(4)+ " is " + norm(W_reweighted(:,5,4),1))
disp("l2 norm of result when epsilon = "+ epsilon_reweighted(4)+ " is " + norm(W_reweighted(:,5,4),2))
%%
% set(gca, 'FontSize', 26);
% xlabel('\theta','fontsize',39);
% ylabel('Intensity(dB)','fontsize',38);
% % title('Beampattern of results with different values of epsilon','fontsize',38);
% % grid on;
% 
% xlim([-18 18])
% ylim([-70 0])
%% Get the readings
close all
clc
clear

% Get the readings, use a larger L
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;
N = 1100;                % to have a aperture of 1020 lambda
d = param.lambda;
L = 200001;             % The resolution is at least 0.035 degree
fov.left = -18;
fov.right = 18;
Res = (fov.right - fov.left)/(L - 1);
angle = linspace(fov.left,fov.right,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
%% LASSO results

load('LASSO_results_and_epsilons_used.mat')
selection_vector = [9 40 59 67];
figure;
a1 = plot(angle,20*log10(abs(S*w_epsilon_increase_all(:,selection_vector(1)))),'LineWidth', 2);
M1 = "\epsilon = " + epsilon_all(selection_vector(1));
% hold on 
figure
a2 = plot(angle,20*log10(abs(S*w_epsilon_increase_all(:,selection_vector(2)))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M1 = "\epsilon = " + epsilon_all(selection_vector(2));
figure
a3 = plot(angle,20*log10(abs(S*w_epsilon_increase_all(:,selection_vector(3)))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M1 = "\epsilon = " + epsilon_all(selection_vector(3));
figure
a4 = plot(angle,20*log10(abs(S*w_epsilon_increase_all(:,selection_vector(4)))),'LineWidth', 2);
M1 = "\epsilon = " + epsilon_all(selection_vector(4));
%%
disp("l1 norm of result when epsilon = "+ epsilon_all(selection_vector(1))+ " is " + norm(w_epsilon_increase_all(:,selection_vector(1)),1))
disp("l2 norm of result when epsilon = "+ epsilon_all(selection_vector(1))+ " is " + norm(w_epsilon_increase_all(:,selection_vector(1)),2))

disp("l1 norm of result when epsilon = "+ epsilon_all(selection_vector(2))+ " is " + norm(w_epsilon_increase_all(:,selection_vector(2)),1))
disp("l2 norm of result when epsilon = "+ epsilon_all(selection_vector(2))+ " is " + norm(w_epsilon_increase_all(:,selection_vector(2)),2))

disp("l1 norm of result when epsilon = "+ epsilon_all(selection_vector(3))+ " is " + norm(w_epsilon_increase_all(:,selection_vector(3)),1))
disp("l2 norm of result when epsilon = "+ epsilon_all(selection_vector(3))+ " is " + norm(w_epsilon_increase_all(:,selection_vector(3)),2))

disp("l1 norm of result when epsilon = "+ epsilon_all(selection_vector(4))+ " is " + norm(w_epsilon_increase_all(:,selection_vector(4)),1))
disp("l2 norm of result when epsilon = "+ epsilon_all(selection_vector(4))+ " is " + norm(w_epsilon_increase_all(:,selection_vector(4)),2))


%% Thresholding results
epsilon_thresholding = epsilon_all(40);
w_2thresholding = w_epsilon_increase_all(:,40);
[W_thresheld, E_thresheld] = threshold(w_2thresholding);

figure
a2 = plot(angle,20*log10(abs(S*W_thresheld(:,601))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "k = " +nnz(W_thresheld(:,601));
figure
a3 = plot(angle,20*log10(abs(S*W_thresheld(:,1001))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M3 = "k = " +nnz(W_thresheld(:,1001));

%% 
disp("l1 norm of result when k = 600 is " + norm(W_thresheld(:,601),1))
disp("l2 norm of result when k = 600 is " + norm(W_thresheld(:,601),2))

disp("l1 norm of result when k = 100 is " + norm(W_thresheld(:,1001),1))
disp("l2 norm of result when k = 100 is " + norm(W_thresheld(:,1001),2))