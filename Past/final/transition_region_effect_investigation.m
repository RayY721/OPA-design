%% transition region investigation; This script shows the effect of using the transition region
% w_modified is the same to w_modified2. w_modified2 is to verify the new
% adjustment function produce the same result as old adjustment function

% The new adjustment function only relax outside the main beam region, in accordance to digital filter design
% The conclusion: the trasition relaxation widen the mainbeam and suppress
% the sidelobe. This term trading off the mainlobe width and the peak
% sidelobe level

% Hoever, the epsilon need to be re-determined after this transition
% relaxation. In this case the epsilon is reduced from 0.4 to 0.06
close all
clc
clear
load('transition_region_effect_investigation.mat')
%% Parameters
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;
N = 1000;                % to have a aperture of 1020 lambda
d = param.lambda;
L = 1000;             % The resolution is at least 0.035 degree
fov.left = -18;
fov.right = 18;
Res = (fov.right - fov.left)/(L - 1);
angle = linspace(fov.left,fov.right,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
Sn = S./sqrt(N);

%% Set the desired beampattern
beam.left = -0.04;
beam.right = 0.04;
desired_pattern = upperboundgen(fov.left,fov.right,[beam.left beam.right],L);
%% Modify the S matrix and reference pattern
loose_range = 0.02;
S_modified = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
desired_pattern_modified = adjustRefpattern2(desired_pattern,[fov.left fov.right],[beam.left beam.right],L,loose_range);

% S_modified2 = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
% desired_pattern_modified2 = adjustRefpattern2(desired_pattern,[fov.left fov.right],[beam.left beam.right],L,loose_range);
%%
cvx_begin
    variable w_modified(N) complex
    minimize(norm(w_modified,1))
    subject to
        norm(S_modified*w_modified - desired_pattern_modified,2) <= 0.06;
cvx_end
%%
% cvx_begin
%     variable w_modified2(N) complex
%     minimize(norm(w_modified2,1))
%     subject to
%         norm(S_modified2*w_modified2 - desired_pattern_modified2,2) <= 0.2;
% cvx_end
%%
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S*w - desired_pattern,2) <= 0.4;
cvx_end
%%
L = 3000;             % The resolution is at least 0.035 degree
fov.left = -18;
fov.right = 18;
Res = (fov.right - fov.left)/(L - 1);
angle = linspace(fov.left,fov.right,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
figure
subplot(2,1,1)
plot(angle,20*log10(abs(S*w_modified)))
ylim([-120 0])
title('The beampattern with transition region used')
ylabel('Intensity(dB)')
xlabel('Angle')
subplot(2,1,2)
plot(angle,20*log10(abs(S*w)))
ylim([-120 0])
title('The beampattern without transition region used')
ylabel('Intensity(dB)')
xlabel('Angle')
figure
subplot(2,1,1)
semilogy(flip(sort(abs(w_modified))))
subplot(2,1,2)
semilogy(flip(sort(abs(w))))