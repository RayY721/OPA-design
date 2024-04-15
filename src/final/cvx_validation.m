%% validation of the cvx algorithm, see if cvx is able to produce the complex solution
% and also why the LS solution is complex.
% starts from a problem with small dimensions

close all
clc
clear

%% Parameters using small value to validate the algorithm
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;
N = 500;                % to have a aperture of 1020 lambda
d = param.lambda;
L = 500;             % The resolution is at least 0.035 degree
fov.left = -18;
fov.right = 18;
Res = (fov.right - fov.left)/(L - 1);
angle = linspace(fov.left,fov.right,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix

%% Set the desired beampattern
beam.left = -3;
beam.right = 3;
p_d = upperboundgen(fov.left,fov.right,[beam.left beam.right],L);
p_d = 1i*p_d;

%% Determine the range of the epsilon
w_ls = pinv(S)*p_d;     % This solution gives the LS error
epsilon_ls = norm(S*w_ls - p_d,2);
epsilon_max = norm(p_d,2);

%%
cvx_begin
    cvx_precision low
    variable w_ls_cvx(N) complex
    minimize(norm(S*w_ls_cvx - p_d,2))
cvx_end

%%  
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S*w - p_d,2) <= 0.6;
cvx_end

%%
figure
plot(angle,20*log10(abs(S*w)))