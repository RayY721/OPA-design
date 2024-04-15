%% This script exam the choice of the epsilon
close all
clc
clear
load('investigation_transition_region_relaxation.mat')
load('investigation_transition_region_relaxation_epsilon_115.mat')

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
%% Modify the S matrix and reference pattern
loose_range = 0.02;
S_modified_002 = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
p_d_modified_002 = adjustRefpattern2(p_d,[fov.left fov.right],[beam.left beam.right],L,loose_range);
% same number of removal as 0.02
% loose_range = 0.03;
% S_modified_003 = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
% p_d_modified_003 = adjustRefpattern2(p_d,[fov.left fov.right],[beam.left beam.right],L,loose_range);
% loose_range = 0.04;
% S_modified_004 = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
% p_d_modified_004 = adjustRefpattern2(p_d,[fov.left fov.right],[beam.left beam.right],L,loose_range);

loose_range = 0.05;
S_modified_005 = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
p_d_modified_005 = adjustRefpattern2(p_d,[fov.left fov.right],[beam.left beam.right],L,loose_range);
% same as 0.05
% loose_range = 0.06;
% S_modified_006 = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
% p_d_modified_006 = adjustRefpattern2(p_d,[fov.left fov.right],[beam.left beam.right],L,loose_range);
% loose_range = 0.07;
% S_modified_007 = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
% p_d_modified_007 = adjustRefpattern2(p_d,[fov.left fov.right],[beam.left beam.right],L,loose_range);

loose_range = 0.08;
S_modified_008 = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
p_d_modified_008 = adjustRefpattern2(p_d,[fov.left fov.right],[beam.left beam.right],L,loose_range);
%% Determine the range of the epsilon
w_ls002 = pinv(S_modified_002)*p_d_modified_002;     % This solution gives the LS error
epsilon_ls002 = norm(S_modified_002*w_ls002 - p_d_modified_002,2);
epsilon_max002 = norm(p_d_modified_002,2);

w_ls005 = pinv(S_modified_005)*p_d_modified_005;     % This solution gives the LS error
epsilon_ls005 = norm(S_modified_005*w_ls005 - p_d_modified_005,2);
epsilon_max005 = norm(p_d_modified_005,2);

w_ls008 = pinv(S_modified_008)*p_d_modified_008;     % This solution gives the LS error
epsilon_ls008 = norm(S_modified_008*w_ls008 - p_d_modified_008,2);
epsilon_max008 = norm(p_d_modified_008,2);

%%
figure

a1 = plot(angle,20*log10(abs(S*w_ls002)),'LineWidth', 2);
M1 = "Transition region = 0.02\circ";
hold on
a2 = plot(angle,20*log10(abs(S*w_ls005)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "Transition region = 0.05\circ";
a3 = plot(angle,20*log10(abs(S*w_ls008)),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M3 = "Transition region = 0.08\circ";
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampattern of LS solutions with different sizes of transition region','fontsize',38);
legend([a1,a2,a3],[M1,M2,M3],'FontSize',39)
xlim([-1,1])
ylim([-90,0])

%% 
cvx_begin
    variable w002(N) complex
    minimize(norm(w002,1))
    subject to
        norm(S_modified_002*w002 - p_d_modified_002,2) <= 0.06;
cvx_end

cvx_begin
    variable w005(N) complex
    minimize(norm(w005,1))
    subject to
        norm(S_modified_005*w005 - p_d_modified_005,2) <= 0.06;
cvx_end

cvx_begin
    variable w008(N) complex
    minimize(norm(w008,1))
    subject to
        norm(S_modified_008*w008 - p_d_modified_008,2) <= 0.06;
cvx_end
%% 
cvx_begin
    variable w002_e115(N) complex
    minimize(norm(w002_e115,1))
    subject to
        norm(S_modified_002*w002_e115 - p_d_modified_002,2) <= 1.15;
cvx_end

cvx_begin
    variable w005_e115(N) complex
    minimize(norm(w005_e115,1))
    subject to
        norm(S_modified_005*w005_e115 - p_d_modified_005,2) <= 1.15;
cvx_end

cvx_begin
    variable w008_e115(N) complex
    minimize(norm(w008_e115,1))
    subject to
        norm(S_modified_008*w008_e115 - p_d_modified_008,2) <= 1.15;
cvx_end

%%
figure;
% plot(angle,20*log10(abs(S*w)))
a1 = plot(angle,20*log10(abs(S*w002)),'LineWidth', 2);
M1 = "Transition region = 0.02\circ";
hold on 
a2 = plot(angle,20*log10(abs(S*w005)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "Transition region = 0.05\circ";
a3 = plot(angle,20*log10(abs(S*w008)),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M3 = "Transition region = 0.08\circ";
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampattern of results with different sizes of transition region with \epsilon = 0.06','fontsize',38);
% grid on;
legend([a1,a2,a3],[M1,M2,M3],'FontSize',39)
xlim([-1 1])
ylim([-60 0])

%% Beampattern with unit power
figure
subplot(1,3,1)
plot(angle,20*log10(abs(S*w002./norm(w002,2))))
xlabel('Angle in degree')
ylabel('Intensity in dB')
% ylim([-60,0])
xlim([-18,18])
title("Transition region = 0.02\circ")
subplot(1,3,2)
plot(angle,20*log10(abs(S*w005./norm(w005,2))))
xlabel('Angle in degree')
ylabel('Intensity in dB')
% ylim([-60,0])
xlim([-18,18])
title("Transition region = 0.02\circ")
subplot(1,3,3)
plot(angle,20*log10(abs(S*w008./norm(w008,2))))
xlabel('Angle in degree')
ylabel('Intensity in dB')
% ylim([-60,0])
xlim([-18,18])
title("Transition region = 0.02\circ")


%%
cvx_begin
    variable w008_001(N) complex
    minimize(norm(w008_001,1))
    subject to
        norm(S_modified_008*w008_001 - p_d_modified_008,2) <= 0.01;
cvx_end
%%
cvx_begin
    variable w005_001(N) complex
    minimize(norm(w005_001,1))
    subject to
        norm(S_modified_005*w005_001 - p_d_modified_005,2) <= 0.01;
cvx_end
%% 
figure
subplot(2,4,1)
plot(angle,20*log10(abs(S*w002)))
xlabel('Angle in degree')
ylabel('Intensity in dB')
% ylim([-60,0])
xlim([-18,18])
title("relaxation = 0.02" )
subplot(2,4,2)
plot(angle,20*log10(abs(S*w005)))
xlabel('Angle in degree')
ylabel('Intensity in dB')
% ylim([-60,0])
xlim([-18,18])
title("relaxation = 0.05" )
subplot(2,4,3)
plot(angle,20*log10(abs(S*w008)))
xlabel('Angle in degree')
ylabel('Intensity in dB')
% ylim([-60,0])
xlim([-18,18])
title("relaxation = 0.08" )
% subplot(2,4,4)
% plot(angle,20*log10(abs(S*w005_0007)))
% xlabel('Angle in degree')
% ylabel('Intensity in dB')
% % ylim([-60,0])
% xlim([-18,18])

subplot(2,4,5)
stem(abs(w002))
xlabel('Index')
ylabel('Amplitude')
% ylim([-60,0])
xlim([0,1100])
title("relaxation = 0.02" )
subplot(2,4,6)
stem(abs(w005))
xlabel('Index')
ylabel('Amplitude')
% ylim([-60,0])
xlim([0,1100])
title("relaxation = 0.05" )
subplot(2,4,7)
stem(abs(w008))
xlabel('Index')
ylabel('Amplitude')
% ylim([-60,0])
xlim([0,1100])
title("relaxation = 0.08" )

% subplot(2,4,8)
% stem(abs(w005_0007))
% xlabel('Index')
% ylabel('Amplitude')
% % ylim([-60,0])
% xlim([0,1100])
% title("relaxation = 0.08 but epsilon smaller" )

%% 
figure
a1 = semilogy(sort(abs(w002),'descend')./norm(w002,2),'LineWidth', 1);
% a1 = semilogy(sort(abs(w_epsilon_increase(:,1)),'descend'),'LineWidth', 1);
M1 = "Transition = 0.02\circ";
hold on
a2 = semilogy(sort(abs(w005),'descend')./norm(w005,2),'Color',[0.9290 0.6940 0.1250],'LineWidth', 1);
% a2 = semilogy(sort(abs(w_epsilon_increase(:,2)),'descend'),'Color',[0.9290 0.6940 0.1250],'LineWidth', 1);
M2 = "Transition = 0.05\circ";
a3 = semilogy(sort(abs(w008),'descend')./norm(w008,2),'Color',[0.4660 0.6740 0.1880],'LineWidth', 1);
% a3 = semilogy(sort(abs(w_epsilon_increase(:,3)),'descend'),'Color',[0.4660 0.6740 0.1880],'LineWidth', 1);
M3 = "Transition = 0.08\circ";
% a4 = semilogy(sort(abs(w_epsilon_increase(:,4)),'descend')./norm(w_epsilon_increase(:,4),2),'LineWidth', 1);
% % a4 = semilogy(sort(abs(w_epsilon_increase(:,4)),'descend'),'LineWidth', 1);
% M4 = "\epsilon = "+epsilon_(4);
set(gca, 'FontSize', 24);
xlabel('The index of elements','fontsize',37);
ylabel('The amplitude of excitation(dB)','fontsize',36);
title('The sorted unit power excitation distribution with different values of epsilon','fontsize',36)
grid on
legend([a1,a2,a3],[M1,M2,M3],'FontSize',39)
xlim([0,1100])