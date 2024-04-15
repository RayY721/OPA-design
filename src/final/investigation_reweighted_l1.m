%% Thresholding investigation
close all
clc
clear
load('investigation_reweighted_l1.mat')

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
%% Determine the range of the epsilon
w_ls = pinv(S)*p_d;     % This solution gives the LS error
epsilon_ls = norm(S*w_ls - p_d,2);
epsilon_max = norm(p_d,2);

%% weighted l1 norm (epsilon = 0.4)
Z = eye(N);
P = 5;
W_rew065 = zeros(N,P);
for i = 1:1:P    
cvx_begin
    variable w_rew(N) complex
    minimize(norm(Z*w_rew,1))
    subject to
        norm(S*w_rew - p_d,2) <= 0.65;       
cvx_end
W_rew065(:,i) = w_rew;
Z = inv(diag(abs(w_rew) + 0.0002));
end

% weighted l1 norm (epsilon = 0.9)
Z = eye(N);
P = 5;
W_rew115 = zeros(N,P);
for i = 1:1:P    
cvx_begin
    variable w_rew(N) complex
    minimize(norm(Z*w_rew,1))
    subject to
        norm(S*w_rew - p_d,2) <= 1.15;       
cvx_end
W_rew115(:,i) = w_rew;
Z = inv(diag(abs(w_rew) + 0.0002));
end

%% Beampattern comparison 
figure
a1 = plot(angle,20*log10(abs(S*W_rew065(:,1))),'LineWidth', 2);
M1 = "No.1 iteration";
hold on 
a2 = plot(angle,20*log10(abs(S*W_rew065(:,5))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M2 = "No.5 iteration";
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampattern of results at different iteration with \epsilon = 0.65','fontsize',38);
% grid on;
legend([a1,a2],[M1,M2],'FontSize',39)
xlim([-18 18])
ylim([-60 0])

%% convergence of results with epsilon = 1.15
figure
a1 = semilogy(sort(abs(W_rew115(:,1))./norm(W_rew115(:,1),2),'descend'),'LineWidth', 2);
M1 = "No.1 iteration";
hold on
a2 = semilogy(sort(abs(W_rew115(:,2))./norm(W_rew115(:,2),2),'descend'),'LineWidth', 2);
M2 = "No.2 iteration";
a3 = semilogy(sort(abs(W_rew115(:,3))./norm(W_rew115(:,3),2),'descend'),'LineWidth', 2);
M3 = "No.3 iteration";
a4 = semilogy(sort(abs(W_rew115(:,4))./norm(W_rew115(:,4),2),'descend'),'LineWidth', 2);
M4 = "No.4 iteration";
a5 = semilogy(sort(abs(W_rew115(:,5))./norm(W_rew115(:,5),2),'descend'),'LineWidth', 2);
M5 = "No.5 iteration";
set(gca, 'FontSize', 26);
xlabel('The index of elements','fontsize',39);
ylabel('The amplitude of excitation(dB)','fontsize',38);
title('The sorted unit power excitation distribution at different iterations (\epsilon = 1.15)','fontsize',38);
legend([a1,a2,a3,a4,a5],[M1,M2,M3,M4,M5],'FontSize',39)
xlim([0,1100])

figure
a1 = stem(abs(W_rew115(:,1)),'LineWidth', 2);
M1 = "No.1 iteration";
hold on
% a2 = stem(abs(W_rew115(:,2))./norm(W_rew115(:,2),2),'LineWidth', 2);
% M2 = "No.2 iteration";
a3 = stem(abs(W_rew115(:,3)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M3 = "No.3 iteration";
% a4 = stem(abs(W_rew115(:,4))./norm(W_rew115(:,4),2),'LineWidth', 2);
% M4 = "No.4 iteration";
a5 = stem(abs(W_rew115(:,5)),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M5 = "No.5 iteration";
grid on
set(gca, 'FontSize', 26);
xlabel('The index of elements','fontsize',39);
ylabel('The amplitude of excitation','fontsize',38);
title('The excitation distribution at different iterations (\epsilon = 1.15) ','fontsize',38);
legend([a1,a3,a5],[M1,M3,M5],'FontSize',39)
xlim([0,1100])

%% convergence of results with epsilon = 0.65
figure
a1 = semilogy(sort(abs(W_rew065(:,1))./norm(W_rew065(:,1),2),'descend'),'LineWidth', 2);
M1 = "No.1 iteration";
hold on
a2 = semilogy(sort(abs(W_rew065(:,2))./norm(W_rew065(:,2),2),'descend'),'LineWidth', 2);
M2 = "No.2 iteration";
a3 = semilogy(sort(abs(W_rew065(:,3))./norm(W_rew065(:,3),2),'descend'),'LineWidth', 2);
M3 = "No.3 iteration";
a4 = semilogy(sort(abs(W_rew065(:,4))./norm(W_rew065(:,4),2),'descend'),'LineWidth', 2);
M4 = "No.4 iteration";
a5 = semilogy(sort(abs(W_rew065(:,5))./norm(W_rew065(:,5),2),'descend'),'LineWidth', 2);
M5 = "No.5 iteration";
set(gca, 'FontSize', 26);
xlabel('The index of elements','fontsize',39);
ylabel('The amplitude of excitation(dB)','fontsize',38);
title('The sorted unit power excitation distribution at different iterations (\epsilon = 0.65)','fontsize',38);
legend([a1,a2,a3,a4,a5],[M1,M2,M3,M4,M5],'FontSize',39)
xlim([0,1100])

figure
a1 = stem(abs(W_rew065(:,1)),'LineWidth', 2);
M1 = "No.1 iteration";
hold on
% a2 = stem(abs(W_rew115(:,2))./norm(W_rew115(:,2),2),'LineWidth', 2);
% M2 = "No.2 iteration";
a3 = stem(abs(W_rew065(:,3)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M3 = "No.3 iteration";
% a4 = stem(abs(W_rew115(:,4))./norm(W_rew115(:,4),2),'LineWidth', 2);
% M4 = "No.4 iteration";
a5 = stem(abs(W_rew065(:,5)),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M5 = "No.5 iteration";
grid on
set(gca, 'FontSize', 26);
xlabel('The index of elements','fontsize',39);
ylabel('The amplitude of excitation','fontsize',38);
title('The excitation distribution at different iterations (\epsilon = 0.65)','fontsize',38);
legend([a1,a3,a5],[M1,M3,M5],'FontSize',39)
xlim([0,1100])

%% Beampattern
figure
a1 = plot(angle,20*log10(abs(S*W_rew115(:,5))),'LineWidth', 2);
M1 = "\epsilon = 1.15";
hold on
a2 = plot(angle,20*log10(abs(S*W_rew065(:,5))),'LineWidth', 2);
M2 = "\epsilon = 0.65";
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampattern of results with different values of epsilon','fontsize',38);
% grid on;
legend([a1,a2],[M1,M2],'FontSize',39)
xlim([-18 18])
ylim([-60 0])

%% weighted l1 norm with transition relaxation (epsilon = 0.06)
Z = eye(N);
P = 5;
W_rew_transi002 = zeros(N,P);
for i = 1:1:P    
cvx_begin
    variable w_rew(N) complex
    minimize(norm(Z*w_rew,1))
    subject to
        norm(S_modified_002*w_rew - p_d_modified_002,2) <= 0.06;       
cvx_end
W_rew_transi002(:,i) = w_rew;
Z = inv(diag(abs(w_rew) + 0.0002));
end
%% weighted l1 norm with transition relaxation (epsilon = 0.08)
Z = eye(N);
P = 5;
W_rew008_transi002 = zeros(N,P);
for i = 1:1:P    
cvx_begin
    variable w_rew(N) complex
    minimize(norm(Z*w_rew,1))
    subject to
        norm(S_modified_002*w_rew - p_d_modified_002,2) <= 0.08;       
cvx_end
W_rew008_transi002(:,i) = w_rew;
Z = inv(diag(abs(w_rew) + 0.0002));
end

%% convergence of results with transition relaxation and epsilon = 0.08
figure
a1 = semilogy(sort(abs(W_rew008_transi002(:,1))./norm(W_rew008_transi002(:,1),2),'descend'),'LineWidth', 2);
M1 = "No.1 iteration";
hold on
a2 = semilogy(sort(abs(W_rew008_transi002(:,2))./norm(W_rew008_transi002(:,2),2),'descend'),'LineWidth', 2);
M2 = "No.2 iteration";
a3 = semilogy(sort(abs(W_rew008_transi002(:,3))./norm(W_rew008_transi002(:,3),2),'descend'),'LineWidth', 2);
M3 = "No.3 iteration";
a4 = semilogy(sort(abs(W_rew008_transi002(:,4))./norm(W_rew008_transi002(:,4),2),'descend'),'LineWidth', 2);
M4 = "No.4 iteration";
a5 = semilogy(sort(abs(W_rew008_transi002(:,5))./norm(W_rew008_transi002(:,5),2),'descend'),'LineWidth', 2);
M5 = "No.5 iteration";
set(gca, 'FontSize', 26);
xlabel('The index of elements','fontsize',39);
ylabel('The amplitude of excitation(dB)','fontsize',38);
title('The sorted unit power excitation distribution at different iterations (\epsilon = 0.08)','fontsize',38);
legend([a1,a2,a3,a4,a5],[M1,M2,M3,M4,M5],'FontSize',39)
xlim([0,1100])

figure
% a1 = stem(abs(W_rew008_transi002(:,1))./norm(W_rew008_transi002(:,1),2),'LineWidth', 2);
a1 = stem(abs(W_rew008_transi002(:,1)),'LineWidth', 2);
M1 = "No.1 iteration";
hold on
% a2 = stem(abs(W_rew115(:,2))./norm(W_rew115(:,2),2),'LineWidth', 2);
% M2 = "No.2 iteration";
a3 = stem(abs(W_rew008_transi002(:,3)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M3 = "No.3 iteration";
% a4 = stem(abs(W_rew115(:,4))./norm(W_rew115(:,4),2),'LineWidth', 2);
% M4 = "No.4 iteration";
a5 = stem(abs(W_rew008_transi002(:,5)),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M5 = "No.5 iteration";
grid on
set(gca, 'FontSize', 26);
xlabel('The index of elements','fontsize',39);
ylabel('The amplitude of excitation','fontsize',38);
title('The excitation distribution at different iterations (\epsilon = 0.08)','fontsize',38);
legend([a1,a3,a5],[M1,M3,M5],'FontSize',39)
xlim([0,1100])

figure
a1 = plot(angle,20*log10(abs(S*W_rew008_transi002(:,1))),'LineWidth', 2);
M1 = "No.1 iteration";
hold on 
a2 = plot(angle,20*log10(abs(S*W_rew008_transi002(:,5))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M2 = "No.5 iteration";
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampattern of results at different iteration with \epsilon = 0.08','fontsize',38);
% grid on;
legend([a1,a2],[M1,M2],'FontSize',39)
xlim([-18 18])
ylim([-60 0])

%% convergence of results with transition relaxation and epsilon = 0.06
figure
a1 = semilogy(sort(abs(W_rew_transi002(:,1))./norm(W_rew_transi002(:,1),2),'descend'),'LineWidth', 2);
M1 = "No.1 iteration";
hold on
a2 = semilogy(sort(abs(W_rew_transi002(:,2))./norm(W_rew_transi002(:,2),2),'descend'),'LineWidth', 2);
M2 = "No.2 iteration";
a3 = semilogy(sort(abs(W_rew_transi002(:,3))./norm(W_rew_transi002(:,3),2),'descend'),'LineWidth', 2);
M3 = "No.3 iteration";
a4 = semilogy(sort(abs(W_rew_transi002(:,4))./norm(W_rew_transi002(:,4),2),'descend'),'LineWidth', 2);
M4 = "No.4 iteration";
a5 = semilogy(sort(abs(W_rew_transi002(:,5))./norm(W_rew_transi002(:,5),2),'descend'),'LineWidth', 2);
M5 = "No.5 iteration";
set(gca, 'FontSize', 26);
xlabel('The index of elements','fontsize',39);
ylabel('The amplitude of excitation(dB)','fontsize',38);
title('The sorted unit power excitation distribution at different iterations (\epsilon = 0.06)','fontsize',38);
legend([a1,a2,a3,a4,a5],[M1,M2,M3,M4,M5],'FontSize',39)
xlim([0,1100])

figure
a1 = stem(abs(W_rew_transi002(:,1)),'LineWidth', 2);
M1 = "No.1 iteration";
hold on
% a2 = stem(abs(W_rew115(:,2))./norm(W_rew115(:,2),2),'LineWidth', 2);
% M2 = "No.2 iteration";
a3 = stem(abs(W_rew_transi002(:,3)),'Color',[0.9290 0.6940 0.1250],'LineWidth', 2);
M3 = "No.3 iteration";
% a4 = stem(abs(W_rew115(:,4))./norm(W_rew115(:,4),2),'LineWidth', 2);
% M4 = "No.4 iteration";
a5 = stem(abs(W_rew_transi002(:,5)),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M5 = "No.5 iteration";
grid on
set(gca, 'FontSize', 26);
xlabel('The index of elements','fontsize',39);
ylabel('The amplitude of excitation','fontsize',38);
title('The excitation distribution at different iterations (\epsilon = 0.06)','fontsize',38);
legend([a1,a3,a5],[M1,M3,M5],'FontSize',39)
xlim([0,1100])

figure
a1 = plot(angle,20*log10(abs(S*W_rew_transi002(:,1))),'LineWidth', 2);
M1 = "No.1 iteration";
hold on 
a2 = plot(angle,20*log10(abs(S*W_rew_transi002(:,5))),'Color',[0.4660 0.6740 0.1880],'LineWidth', 2);
M2 = "No.5 iteration";
set(gca, 'FontSize', 26);
xlabel('\theta','fontsize',39);
ylabel('Intensity(dB)','fontsize',38);
title('Beampattern of results at different iteration with \epsilon = 0.06','fontsize',38);
% grid on;
legend([a1,a2],[M1,M2],'FontSize',39)
xlim([-18 18])
ylim([-60 0])


%% thresholding 
[W,E] = threshold(W_rew_transi002(:,5));       % subsitude the w with others

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

