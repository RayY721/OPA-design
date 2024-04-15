%% Improve the algorithm a bit
close all
clc
clear

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
steering_vector_0 = exp(1i*param.k*sin(0)*x');      

%% Set the desired beampattern
beam.left = -0.035;
beam.right = 0.035;
p_d = upperboundgen(fov.left,fov.right,[beam.left beam.right],L);

%% CVX for LS solution

cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S* - p_d,2) <= 1;
        steering_vector_0*w == 1;

cvx_end


%% 

plot(20*log10(abs(S*w)))





