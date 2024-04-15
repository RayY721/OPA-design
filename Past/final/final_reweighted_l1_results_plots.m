%% Reweighted l1-norm, with corresponding epsilon values
close all
clc
clear
load('LASSO_results_and_epsilons_used.mat')
selection_vector = [9 40 59 67];
epsilon_reweighted = epsilon_all(selection_vector);
load('reweighted_results.mat')

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

%% weighted l1 norm (epsilon = 0.4)
% Use a tensor to store the results of reweighted algo, 
% For single epsilon value, the result is a (N,P) matrix 
% For J = size(epsilon_reweighted,2) number of epsilon values, the result
% should be an order 3 tensor: (N,P,J) 
J = size(epsilon_reweighted,2);
P = 5;
W_reweighted = zeros(N,P,J);
optimal_value = zeros(P,J);

for j =1:1:J
    Z = eye(N);
    for i = 1:1:P                           % The 
    cvx_begin
        variable w(N) complex
        minimize(norm(Z*w,1))
        subject to
            norm(S*w - p_d,2) <= epsilon_reweighted(j);       
    cvx_end
    W_reweighted(:,i,j) = w;
    Z = inv(diag(abs(w) + 0.0002));
    optimal_value(i,j) = cvx_optval;
    end
end

%% Proof of convergence (the value of objective function) 

figure
a1 = plot(1:5,optimal_value(:,1), '-+','LineWidth',2,'MarkerSize',8);
M1 = "\epsilon = "+epsilon_reweighted(1);
hold on
a2 = plot(1:5,optimal_value(:,2), '-diamond','Color',[0.9290 0.6940 0.1250],'LineWidth',2,'MarkerSize',8);
M2 = "\epsilon = "+epsilon_reweighted(2);
a3 = plot(1:5,optimal_value(:,3), '-square','Color',[0.4660 0.6740 0.1880],'LineWidth',2,'MarkerSize',8);
M3 = "\epsilon = "+epsilon_reweighted(3);
a4 = plot(1:5,optimal_value(:,4), '-x','LineWidth',2,'MarkerSize',8);
M4 = "\epsilon = "+epsilon_reweighted(4);
xticks(1:5)
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
set(gca, 'FontSize', 24);
xlabel('Iteration (m)' ,'fontsize',37);
ylabel('$\|\mathbf{\Gamma}^{(m)}\textbf{\textit{w}}\|_1$','Interpreter','latex','fontsize',36);
title('The Change of Optimal Value of Objective Function')

%% Proof of convergence
for aaa = 1
figure
a1 = semilogy(sort(abs(W_reweighted(:,1,1)),'descend'),'LineWidth', 2);
M1 = "No.1 iteration";
hold on
a2 = semilogy(sort(abs(W_reweighted(:,2,1)),'descend'),'LineWidth', 2);
M2 = "No.2 iteration";
a3 = semilogy(sort(abs(W_reweighted(:,3,1)),'descend'),'LineWidth', 2);
M3 = "No.3 iteration";
a4 = semilogy(sort(abs(W_reweighted(:,4,1)),'descend'),'LineWidth', 2);
M4 = "No.4 iteration";
a5 = semilogy(sort(abs(W_reweighted(:,5,1)),'descend'),'LineWidth', 2);
M5 = "No.5 iteration";
set(gca, 'FontSize', 26);
xlabel('The index of elements','fontsize',39);
ylabel('The amplitude of excitation(dB)','fontsize',38);
title("The sorted excitation distribution at different iterations \epsilon" + epsilon_reweighted(1),'fontsize',38)
legend([a1,a2,a3,a4,a5],[M1,M2,M3,M4,M5],'FontSize',39)
xlim([0,1100])

%%
figure
a1 = semilogy(sort(abs(W_reweighted(:,1,2)),'descend'),'LineWidth', 2);
M1 = "No.1 iteration";
hold on
a2 = semilogy(sort(abs(W_reweighted(:,2,2)),'descend'),'LineWidth', 2);
M2 = "No.2 iteration";
a3 = semilogy(sort(abs(W_reweighted(:,3,2)),'descend'),'LineWidth', 2);
M3 = "No.3 iteration";
a4 = semilogy(sort(abs(W_reweighted(:,4,2)),'descend'),'LineWidth', 2);
M4 = "No.4 iteration";
a5 = semilogy(sort(abs(W_reweighted(:,5,2)),'descend'),'LineWidth', 2);
M5 = "No.5 iteration";
set(gca, 'FontSize', 26);
xlabel('The index of elements','fontsize',39);
ylabel('Amplitude(dB)','fontsize',38);
title("The Sorted Amplitude Distribution of Excitations at different iterations \epsilon = " + epsilon_reweighted(2),'fontsize',38)
legend([a1,a2,a3,a4,a5],[M1,M2,M3,M4,M5],'FontSize',39)
xlim([0,1100])
%%
figure
a1 = semilogy(sort(abs(W_reweighted(:,1,3)),'descend'),'LineWidth', 2);
M1 = "No.1 iteration";
hold on
a2 = semilogy(sort(abs(W_reweighted(:,2,3)),'descend'),'LineWidth', 2);
M2 = "No.2 iteration";
a3 = semilogy(sort(abs(W_reweighted(:,3,3)),'descend'),'LineWidth', 2);
M3 = "No.3 iteration";
a4 = semilogy(sort(abs(W_reweighted(:,4,3)),'descend'),'LineWidth', 2);
M4 = "No.4 iteration";
a5 = semilogy(sort(abs(W_reweighted(:,5,3)),'descend'),'LineWidth', 2);
M5 = "No.5 iteration";
set(gca, 'FontSize', 26);
xlabel('The index of elements','fontsize',39);
ylabel('The amplitude of excitation(dB)','fontsize',38);
title("The sorted excitation distribution at different iterations \epsilon" + epsilon_reweighted(3),'fontsize',38)
legend([a1,a2,a3,a4,a5],[M1,M2,M3,M4,M5],'FontSize',39)
xlim([0,1100])
%%
figure
a1 = semilogy(sort(abs(W_reweighted(:,1,4)),'descend'),'LineWidth', 2);
M1 = "No.1 iteration";
hold on
a2 = semilogy(sort(abs(W_reweighted(:,2,4)),'descend'),'LineWidth', 2);
M2 = "No.2 iteration";
a3 = semilogy(sort(abs(W_reweighted(:,3,4)),'descend'),'LineWidth', 2);
M3 = "No.3 iteration";
a4 = semilogy(sort(abs(W_reweighted(:,4,4)),'descend'),'LineWidth', 2);
M4 = "No.4 iteration";
a5 = semilogy(sort(abs(W_reweighted(:,5,4)),'descend'),'LineWidth', 2);
M5 = "No.5 iteration";
set(gca, 'FontSize', 26);
xlabel('The index of elements','fontsize',39);
ylabel('The amplitude of excitation(dB)','fontsize',38);
title("The sorted excitation distribution at different iterations \epsilon" + epsilon_reweighted(4),'fontsize',38)
legend([a1,a2,a3,a4,a5],[M1,M2,M3,M4,M5],'FontSize',39)
xlim([0,1100])
end
%% Plot the spatial excitation distribution

figure
a1 = stem(abs(W_reweighted(:,1,2)),'LineWidth', 2);
M1 = "k = " +nnz(W_reweighted(:,1));
hold on 
% a2 = stem(abs(W_reweighted(:,2,2)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
% M2 = "k = " +nnz(W_reweighted(:,601));
a3 = stem(abs(W_reweighted(:,5,2)),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M3 = "k = " +nnz(W_reweighted(:,1001));
set(gca, 'FontSize', 24);
xlabel('Index','fontsize',37);
ylabel('|\mathbf{w}_k|','fontsize',36);
title('The excitation distribution with different values of epsilon','fontsize',36)
% legend([a1,a3],[M1,M3],'FontSize',39)
xlim([0,1100])

%% Plot the beampattern

figure
a1 = plot(angle,20*log10(abs(S*W_reweighted(:,1,2))),'LineWidth', 2);
M1 = "No.1 iteration";
hold on
a2 = plot(angle,20*log10(abs(S*W_reweighted(:,5,2))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "No.5 iteration";
% a3 = plot(angle,20*log10(abs(S*W_reweighted(:,5,2))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
% M3 = "k = " ;
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampatterns at the First and the Last Iteration','fontsize',38);
% grid on;
legend([a1,a2],[M1,M2],'FontSize',39)
xlim([-18 18])
ylim([-70 0])
%% Plot the beampattern seperately

figure
plot(angle,20*log10(abs(S*W_reweighted(:,1,2))),'LineWidth', 2);
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampatterns of LASSO','fontsize',38);
% grid on;
xlim([-18 18])
ylim([-70 0])



figure
plot(angle,20*log10(abs(S*W_reweighted(:,5,2))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);

set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampatterns of reweighted l1 norm minimization','fontsize',38);
% grid on;
xlim([-18 18])
ylim([-70 0])


%% Plot the sorted amplitude distribution for different epsilon
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


%% plot the spatial excitation distribution
figure
a1 = stem(abs(W_reweighted(:,5,1)),'LineWidth', 1);
M1 = "\epsilon = " + epsilon_reweighted(1);
hold on 
a2 = stem(abs(W_reweighted(:,5,2)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 1);
M2 = "\epsilon = " + epsilon_reweighted(2);
a3 = stem(abs(W_reweighted(:,5,3)),'Color',[0.4660 0.6740 0.1880],'LineWidth', 1);
M3 = "\epsilon = " + epsilon_reweighted(3);
a4 = stem(abs(W_reweighted(:,5,4)),'LineWidth', 1);
M4 = "\epsilon = " + epsilon_reweighted(4);
set(gca, 'FontSize', 24);
xlabel('Index','fontsize',37);
ylabel('Amplitude','fontsize',36);
% title('The excitation distribution with different values of epsilon','fontsize',36)
legend([a3,a4],[M3,M4],'FontSize',39)
xlim([0,1100])
ylim([0,4e-3])

%% plot the spatial excitation distribution in seperate plot
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

%% The beampattern in one 
figure;
% plot(angle,20*log10(abs(S*w)))
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
% grid on;
legend([a1,a2,a3,a4],[M1,M2,M3,M4],'FontSize',39)
xlim([-18 18])
ylim([-70 0])

%% The investigation of how epsilon is used on the beampattern matching
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







