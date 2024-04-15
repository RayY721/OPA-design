%% Thresholding investigation
close all
clc
clear
load('investigation_transition_region_relaxation.mat')
load('epsilon_investigation_with_range.mat')

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

%% thresholding 
[W_LASSO,E_LASSO] = threshold(w_epsilon_increase(:,1));       % subsitude the w with others
[W_transition,E_transition] = threshold(w008);
%% Animation to observe the effects of thresholding
fig = figure();
% Define the range for the animation
numFrames = 1000;

% Loop through each frame and update the plot
for i = 1:5:numFrames
    % Clear the previous plot
    clf
    % Update the plot
    subplot(2,1,1)
    plot(angle,20*log10(abs(S*W(:,i))))
    % yline(10*log10(MSE_app(i)))
    xlim([-18, 18]);
    ylim([-60, 0]);
    title("Beampattern with k =" + nnz(W(:,i)));
    % xlabel('X-axis');
    % ylabel('Y-axis');

    subplot(2,1,2)
    stem(abs(W(:,i)))
    xlim([0,N])

    % subplot(4,1,2)
    % plot(flip(MSE_new))
    % set(gca,'xdir','reverse')
    % xline(nnze(i))
    % title('MSE exact value')
    % text(nnze(i),MSE_new(i),"\leftarrow MSE = " +MSE_new(i))
    % 
    % subplot(4,1,3)
    % plot(flip(MSE_app))
    % set(gca,'xdir','reverse')
    % xline(nnze(i))
    % title('MSE approximation')
    % text(nnze(i),MSE_app(i),"\leftarrow MSE = " +MSE_app(i))

    % subplot(2,1,2)
    % plot(k,flip(nnze));
    % set(gca,'xdir','reverse')
    % xline(nnze(i))
    
    % Pause for a short time to control animation speed
    pause(0.01);
    % Capture the frame for the animation
    frame = getframe(fig);
    % Convert the frame to an image
    image = frame2im(frame);
    % Write the image to a video file (requires VideoWriter)
    if i == 1
        videoFileName = 'animation_for_w_1.avi';
        vidObj = VideoWriter(videoFileName);
        open(vidObj);
    end
    writeVideo(vidObj, image);
end
% Close the video file
close(vidObj);
disp('Animation created successfully.');


%% compare w_epsilon_increase(:,1) with different sparsity level

figure
a1 = plot(angle,20*log10(abs(S*W_LASSO(:,1))),'LineWidth', 2);
M1 = "k = " +nnz(W_LASSO(:,1));
hold on
a2 = plot(angle,20*log10(abs(S*W_LASSO(:,551))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "k = " +nnz(W_LASSO(:,551));
a3 = plot(angle,20*log10(abs(S*W_LASSO(:,721))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M3 = "k = " +nnz(W_LASSO(:,721));

set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampattern of thresholding solution at different sparsity level k','fontsize',38);
% grid on;
legend([a1,a2,a3],[M1,M2,M3],'FontSize',39)
xlim([-18 18])
ylim([-70 0])

%% compare w_epsilon_increase(:,1) with different sparsity level

figure
a1 = plot(angle,20*log10(abs(S*W_transition(:,1))),'LineWidth', 2);
M1 = "k = " +nnz(W_transition(:,1));
hold on
a2 = plot(angle,20*log10(abs(S*W_transition(:,551))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "k = " +nnz(W_transition(:,551));
a3 = plot(angle,20*log10(abs(S*W_transition(:,721))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M3 = "k = " +nnz(W_transition(:,721));

set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampattern of thresholding solution at different sparsity level k','fontsize',38);
% grid on;
legend([a1,a2,a3],[M1,M2,M3],'FontSize',39)
xlim([-18 18])
ylim([-70 0])

%% compare beampattern w008 and w_epsilon_increase
figure
a1 = plot(angle,20*log10(abs(S*W_LASSO(:,721))),'LineWidth', 2);
M1 = "LASSO solution";
hold on
a2 = plot(angle,20*log10(abs(S*W_transition(:,721))),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M2 = "Modified LASSO solution";
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampattern of thresholding different solutions at same sparsity level k = 380','fontsize',38);
% grid on;
legend([a1,a2],[M1,M2],'FontSize',39)
xlim([-18 18])
ylim([-70 0])

%% compare the spatical distribution of excitation 
figure
a1 = stem(abs(W_LASSO(:,721)),'LineWidth', 1);
M1 = "LASSO solution";
hold on 
a2 = stem(abs(W_transition(:,721)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 1);
M2 = "Modified LASSO solution";
set(gca, 'FontSize', 24);
xlabel('The index of elements','fontsize',37);
ylabel('The amplitude of excitation','fontsize',36);
title('The excitation distribution of thresholding different solutions at same sparsity level k =380','fontsize',36)
legend([a1,a2],[M1,M2],'FontSize',39)
xlim([0,1100])
