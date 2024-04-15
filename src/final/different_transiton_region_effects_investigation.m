%% transition region effects with different values
% conclusion: the larger the transition region is, the lower the peak
% sidelobe is and wider the mainbeam width is. 
% However, the epsilon requires further tightening as well. Simply relaxing
% the trasition region without tightening epsilon makes the solution less sparse
close all
clc
clear
load('different_transiton_region_effects_investigation.mat')
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
% Set the desired beampattern
beam.left = -0.04;
beam.right = 0.04;
desired_pattern = upperboundgen(fov.left,fov.right,[beam.left beam.right],L);
%% Modify the S matrix and reference pattern
loose_range = 0.02;
S_modified = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
desired_pattern_modified = adjustRefpattern2(desired_pattern,[fov.left fov.right],[beam.left beam.right],L,loose_range);

cvx_begin
    variable w_loose2(N) complex
    minimize(norm(w_loose2,1))
    subject to
        norm(S_modified*w_loose2 - desired_pattern_modified,2) <= 0.06;
cvx_end
%%
loose_range = 0.06;
S_modified = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
desired_pattern_modified = adjustRefpattern2(desired_pattern,[fov.left fov.right],[beam.left beam.right],L,loose_range);

cvx_begin
    variable w_loose6(N) complex
    minimize(norm(w_loose6,1))
    subject to
        norm(S_modified*w_loose6 - desired_pattern_modified,2) <= 0.01;
cvx_end
%%
loose_range = 0.09;
S_modified = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
desired_pattern_modified = adjustRefpattern2(desired_pattern,[fov.left fov.right],[beam.left beam.right],L,loose_range);

cvx_begin
    variable w_loose9(N) complex
    minimize(norm(w_loose9,1))
    subject to
        norm(S_modified*w_loose9 - desired_pattern_modified,2) <= 0.01;
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
subplot(3,1,1)
plot(angle,20*log10(abs(S*w_loose2)))
ylim([-120 0])
title('The beampattern with transition region of 0.02 degree')
ylabel('Intensity(dB)')
xlabel('Angle')
subplot(3,1,2)
plot(angle,20*log10(abs(S*w_loose6)))
ylim([-120 0])
title('The beampattern with transition region of 0.06 degree')
ylabel('Intensity(dB)')
xlabel('Angle')
subplot(3,1,3)
plot(angle,20*log10(abs(S*w_loose9)))
ylim([-120 0])
title('The beampattern with transition region of 0.09 degree')
ylabel('Intensity(dB)')
xlabel('Angle')