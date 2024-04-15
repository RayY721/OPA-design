%% script for mid-term
close all
clc
clear
%% Parameters
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;
N = 1000;                % to have a aperture of 1020 lambda
d = param.lambda;
L = 1000;             % now Lz = 700 or 900 doesn't work
win.leftend = -18;
win.rightend = 18;
Res = (win.rightend - win.leftend)/(L - 1);
angle = linspace(win.leftend,win.rightend,L)';  % in degree
theta = angle*pi/180;
angle90 = linspace(-90,90,3*L)';
theta90 = angle90*pi/180;
% u = linspace(sin(win.leftend*pi/180),sin(win.rightend*pi/180),L)';
x = (-(N-1)*d/2:d:(N-1)*d/2)';
% x = (0:d:(N-1)*d)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
S90 = exp(1i*param.k*sin(theta90)*x');  
% S_u = exp(1i*param.k*u*x');
%% Set the desired beam pattern (Loose than our requirement) 
% left_win = -0.4;
% right_win = 0.4;
left_win = -0.04;
right_win = 0.04;

desired_pattern = upperboundgen(win.leftend,win.rightend,[left_win right_win],L);       % The desired pattern
figure
plot(angle,desired_pattern)
xlabel('Angle (degrees)');
ylabel('Amplitude')
title('Ideal beam pattern')
%% Modify the S matrix and reference pattern
% loose_range = 0.2;    % in degree   (0.014 is not feasible)
loose_range = 0.02;

desired_pattern_modified = adjustrefpattern(desired_pattern,[win.leftend win.rightend],[left_win right_win],L,loose_range);

S_new = adjustSmatrix(S,[win.leftend win.rightend],[left_win right_win],L,loose_range);
%% Instead of runnning the following sections, load the result directly. The result is stored at the "W_4_midterm.mat"
load("W_4_midterm.mat")
P = 5;
%% cvx 
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S_new*w - desired_pattern_modified,2) <= 0.1;       % 0.1 for L=1000, 0.72 for L=10000
        % norm(S_newn*wn - desired_pattern_with_noise_modified,2) <= 0.1;      
cvx_end
%% weighted l1 norm (epsilon = 0.1)
Z = eye(N);
P = 5;
W_rew = zeros(N,P);
for i = 1:1:P    
cvx_begin
    variable w_rew(N) complex
    minimize(norm(Z*w_rew,1))
    subject to
        norm(S_new*w_rew - desired_pattern_modified,2) <= 0.1;       % 0.1 for L=1000, 0.72 for L=10000   
cvx_end
W_rew(:,i) = w_rew;
Z = inv(diag(abs(w_rew) + 0.0002));
end
%% weighted l1 norm (epsilon = 0.2
Z = eye(N);
P = 5;
W_rew2 = zeros(N,P);
for i = 1:1:P    
cvx_begin
    variable w_rew2(N) complex
    minimize(norm(Z*w_rew2,1))
    subject to
        norm(S_new*w_rew2 - desired_pattern_modified,2) <= 0.2;       % 0.1 for L=1000, 0.72 for L=10000     
cvx_end
W_rew2(:,i) = w_rew2;
Z = inv(diag(abs(w_rew2) + 0.0002));
end
%% The following sections should always be run 
semilogy(flip(sort(abs(w))))
title('Sorted amplitude distribution of w')
ylabel('Amplitude(dB)')
xlabel('Index')
figure
plot(angle,20*log10(abs(S*w)))
title('Array pattern')
ylabel('Intensity(dB)')
xlabel('Angle')
figure
plot(angle90,20*log10(abs(S90*w)))
title('Array pattern')
ylabel('Intensity(dB)')
xlabel('Angle')
figure
stem(abs(w));
title('Amplitude distribution')
ylabel('Amplitude')
xlabel('Index')
%% plot the weighted result
figure;
for i = 1:1:P
subplot(P,1,i)
semilogy(flip(sort(abs(W_rew(:,i)))))
end
%%
figure;
for i = 1:1:P
subplot(P,1,i)
plot(angle,20*log10(abs(S*W_rew(:,i))))
end
%% compare result between l1 norm and weighted l1 norm
figure
subplot(2,1,1)
semilogy(flip(sort(abs(W_rew(:,1)))))
title('Sorted distribution of w')
ylabel('Amplitude')
xlabel('Index')
subplot(2,1,2)
semilogy(flip(sort(abs(W_rew2(:,5)))))
title('Sorted distribution of w')
ylabel('Amplitude')
xlabel('Index')
%% compare pattern between l1 norm and weighted l1 norm
figure
subplot(2,1,1)
plot(angle,20*log10(abs(S*W_rew(:,1))))
title('Original result pattern')
xlabel('Angle (degrees)');
ylabel('Intensity(dB)')
subplot(2,1,2)
plot(angle,20*log10(abs(S*W_rew(:,5))))
title('Enhanced sparsity pattern')
xlabel('Angle (degrees)');
ylabel('Intensity(dB)')

%% plot the weighted result
figure;
for i = 1:1:P
subplot(P,1,i)
semilogy(flip(sort(abs(W_rew2(:,i)))))
end
%%
figure;
for i = 1:1:P
subplot(P,1,i)
plot(angle,20*log10(abs(S*W_rew2(:,i))))
end
%% compare result and pattern between weighted l1 norm with different epsilon
figure
subplot(2,1,1)
semilogy(flip(sort(abs(W_rew(:,5)))))
title('Sorted distribution of w')
ylabel('Amplitude')
xlabel('Index')
subplot(2,1,2)  
semilogy(flip(sort(abs(W_rew2(:,5)))))
title('Sorted distribution of w')
ylabel('Amplitude')
xlabel('Index')

figure
subplot(2,1,1)
plot(angle,20*log10(abs(S*W_rew(:,5))))
title('Original result pattern')
xlabel('Angle (degrees)');
ylabel('Intensity(dB)')
subplot(2,1,2)
plot(angle,20*log10(abs(S*W_rew2(:,5))))
title('Enhanced sparsity pattern')
xlabel('Angle (degrees)');
ylabel('Intensity(dB)')
%% Thresholding the original result
% [cost,th_level,w_th,W,MSE,MSE_new,MSE_app] = iterativethreshold(w,S,0.002);
[W,E] = threshold(w);
%% find best result with a given peak sidelobe level
W_result = findbestsolution(W,S,-25,[win.leftend win.rightend],[-1 1]);
figure
plot(angle,20*log10(abs(S*W_result)))
xline(-1)
xline(1)
title("Array pattern with k =" + nnz(W_result))
%% Animation to show that the pattern and error in w and error in pattern
% Create a figure
fig = figure();
% Define the range for the animation
numFrames = 1000;
k_vector = 1000:-1:1;
sortedw = flip(sort(abs(w)));
% Loop through each frame and update the plot
for i = 1:20:numFrames
    % Clear the previous plot
    clf
    % Update the plot
    subplot(3,1,1)
    plot(angle,20*log10(abs(S*W(:,N-i))))
    % yline(10*log10(MSE_app(i)))
    % xline(left_win)
    % xline(right_win)
    xlim([-18, 18]);
    ylim([-80, 0]);
    title("Array pattern with k =" + nnz(W(:,N-i)));
    xlabel('Angle (degrees)');
    ylabel('Intensity(dB)')
    % subplot(3,1,2)
    % plot(abs(E(:,N-i)));
    % ylim([0, 6e-3])
    % title("Difference introduced in w")
    % xlabel("Index")
    % % set(gca,'xdir','reverse')
    % % xline(nnze(i))
    subplot(3,1,2)
    plot(angle,20*log10(abs(S*E(:,N-i))))
    ylim([-100,0])
    title("Difference introduced on pattern")
    xlabel('Angle (degrees)');
    xlim([-18, 18]);
    ylabel('Intensity(dB)')
    
    subplot(3,1,3)
    semilogy(sortedw)
    title("Sorted distribution of w")
    text(k_vector(i),sortedw(N-i),'\leftarrow threshold')
    xline(k_vector(i))
    xlabel('Index')
    ylabel('Amplitude(dB)')

    % Pause for a short time to control animation speed
    pause(0.01);
    % Capture the frame for the animation
    frame = getframe(fig);
    % Convert the frame to an image
    image = frame2im(frame);
    % Write the image to a video file (requires VideoWriter)
    if i == 1
        videoFileName = 'change_of_pattern.mp4';
        vidObj = VideoWriter(videoFileName,'MPEG-4');
        vidObj.FrameRate = 10;
        open(vidObj);
    end
    writeVideo(vidObj, image);
end

% Close the video file
close(vidObj);

disp('Animation created successfully.');

%% Thresholding the enhanced result
[Wrew2,Erew2] = threshold(W_rew2(:,5));

%% find best result with a given peak sidelobe level
W_result_enchanced = findbestsolution(Wrew2,S,-25,[win.leftend win.rightend],[-1 1]);
figure
plot(angle,20*log10(abs(S*W_result_enchanced)))
% xline(-1)
% xline(1)
title("Array pattern with k =" + nnz(W_result_enchanced))
%% Animation to show that the pattern and error in enhanced w and error in pattern
% Create a figure
fig = figure();

% Define the range for the animation
numFrames = 1000;
theta = linspace(0, 2*pi, numFrames);
k_vector = 1000:-1:1;
sortedw_rew2 = flip(sort(abs(W_rew2(:,5))));
% Loop through each frame and update the plot
for i = 500:11:numFrames
    % Clear the previous plot
    clf
    % Update the plot
    subplot(3,1,1)
    plot(angle,20*log10(abs(S*Wrew2(:,N-i))))
    % yline(10*log10(MSE_app(i)))
    % xline(left_win)
    % xline(right_win)
    xlim([-18, 18]);
    ylim([-80, 0]);
    title("Array pattern with k =" + nnz(Wrew2(:,N-i)));
    xlabel('Angle (degrees)');
    ylabel('Intensity(dB)')
    % subplot(4,1,2)
    % plot(abs(Erew2(:,N-i)));
    % ylim([0, 6e-3])
    % title("Difference introduced in w")
    % xlabel("Index")
    % set(gca,'xdir','reverse')
    % xline(nnze(i))
    subplot(3,1,2)
    plot(angle,20*log10(abs(S*Erew2(:,N-i))))
    ylim([-100,0])
    title("Difference introduced on pattern")
    xlabel('Angle (degrees)');
    xlim([-18, 18]);
    ylabel('Intensity(dB)')
    
    subplot(3,1,3)
    semilogy(sortedw_rew2)
    xline(k_vector(i))
    text(k_vector(i),sortedw_rew2(N-i),'\leftarrow threshold')
    xlabel('Index')
    ylabel('Amplitude(dB)')
    title("Sorted distribution of w")

    % Pause for a short time to control animation speed
    pause(0.01);
    % Capture the frame for the animation
    frame = getframe(fig);
    % Convert the frame to an image
    image = frame2im(frame);
    % Write the image to a video file (requires VideoWriter)
    if i == 500
        videoFileName = 'change_of_pattern_enchanced2.mp4';
        vidObj = VideoWriter(videoFileName,'MPEG-4');
        vidObj.FrameRate = 10;
        open(vidObj);
    end
    writeVideo(vidObj, image);
end

% Close the video file
close(vidObj);

disp('Animation created successfully.');

%% Display the spacial distribution of active sensor
nzindex = find(W_result_enchanced);
plot(x(nzindex)/param.lambda,zeros(length(nzindex),1),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8)
grid on
%%
nzindex_shifted = zeros(size(nzindex));
nzindex_shifted(2:end) = nzindex(1:end-1);
spacing_in_lambda = nzindex(2:end) - nzindex_shifted(2:end);
figure
stem(spacing_in_lambda)
title('The spacing between previous element (counting from the left)')
xlabel('Index')
ylabel('Spacing in wavelength')
figure
histogram(spacing_in_lambda)
title('Histogram of the spacing')
xlabel('The number of spacing from the first active element at the left side')
ylabel('Spacing between active element in d')
%%
figure
stem(abs(W_result_enchanced))
xlabel('Index of element')
ylabel('Amplitude of element')
title('Amplitude distribution of result')
%% Display the phase distribution of active sensor
stem(angle(W_result_enchanced));
% = zeros(size(nzindex));
% Left it here for now due to the limited time. But the idea is that to
% inspect the solution from its amplitude distribution and phase
% distribution. 