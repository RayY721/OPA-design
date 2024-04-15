%% This script investigate the reweighted l1 norm
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

%% weighted l1 norm (epsilon = 0.06)
Z = eye(N);
P = 5;
W_rew = zeros(N,P);
for i = 1:1:P    
cvx_begin
    variable w_rew(N) complex
    minimize(norm(Z*w_rew,1))
    subject to
        norm(S_modified*w_rew - desired_pattern_modified,2) <= 0.06;       
cvx_end
W_rew(:,i) = w_rew;
Z = inv(diag(abs(w_rew) + 0.0002));
end

%% plot the weighted result
figure;
for i = 1:1:P
subplot(P,1,i)
semilogy(flip(sort(abs(W_rew(:,i)))))
xlabel('Index')
ylabel('Amplitude(dB)')
end
sgtitle('The sorted amplitude distribution of excitations with epsilon = 0.06')
%%
figure;
for i = 1:1:P
subplot(P,1,i)
plot(angle,20*log10(abs(S*W_rew(:,i))))
end
%% weighted l1 norm (epsilon = 0.09)
Z = eye(N);
P = 5;
W_rew2 = zeros(N,P);
for i = 1:1:P    
cvx_begin
    variable w_rew2(N) complex
    minimize(norm(Z*w_rew2,1))
    subject to
        norm(S_modified*w_rew2 - desired_pattern_modified,2) <= 0.09;       
cvx_end
W_rew2(:,i) = w_rew2;
Z = inv(diag(abs(w_rew2) + 0.0002));
end
%% plot the weighted result
figure;
for i = 1:1:P
subplot(P,1,i)
semilogy(flip(sort(abs(W_rew2(:,i)))))
xlabel('Index')
ylabel('Amplitude(dB)')
end
sgtitle('The sorted amplitude distribution of excitations with epsilon = 0.09')
%%
figure;
for i = 1:1:P
subplot(P,1,i)
plot(angle,20*log10(abs(S*W_rew2(:,i))))
end
%%
figure
subplot(2,1,1)
semilogy(flip(sort(abs(W_rew(:,1)))))
title('The sorted amplitude distribution of excitations')
xlabel('Index')
ylabel('Amplitude(dB)')
subplot(2,1,2)
semilogy(flip(sort(abs(W_rew(:,end))))) 
title('The sorted amplitude distribution of excitations')
xlabel('Index')
ylabel('Amplitude(dB)')
%%
L = 5000;             % The resolution is at least 0.035 degree
fov.left = -18;
fov.right = 18;
Res = (fov.right - fov.left)/(L - 1);
angle = linspace(fov.left,fov.right,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
figure
subplot(1,2,1)
plot(angle,20*log10(abs(S*W_rew(:,1))))
ylim([-120 0])
title('The beampattern of result of l_1 norm solution')
ylabel('Intensity(dB)')
xlabel('Angle')
subplot(1,2,2)
plot(angle,20*log10(abs(S*W_rew(:,end))))
ylim([-120 0])
title('The beampattern of result of reweighted l_1 norm solution')
ylabel('Intensity(dB)')
xlabel('Angle')