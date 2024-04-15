%% Investigate transition relaxation in LASSO
close all
clc
clear

%% Parameters
for j = 1:1
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
end

%% Set the transition relaxation
loose_range = 0.02;
S_002 = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
p_d_002 = adjustRefpattern2(p_d,[fov.left fov.right],[beam.left beam.right],L,loose_range);

loose_range = 0.05;
S_005 = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
p_d_005 = adjustRefpattern2(p_d,[fov.left fov.right],[beam.left beam.right],L,loose_range);

loose_range = 0.08;
S_008 = adjustSmatrix2(S,[fov.left fov.right],[beam.left beam.right],L,loose_range);
p_d_008 = adjustRefpattern2(p_d,[fov.left fov.right],[beam.left beam.right],L,loose_range);

%% LS solution with pinv()

w_ls = pinv(S)*1000*p_d;                 % This solution gives the LS error
epsilon_ls = norm(S*w_ls/1000 - p_d,2);
epsilon_max = norm(p_d,2);          % The upper bound of the 
%%
w_ls = pinv(S)*p_d;                 % This solution gives the LS error
epsilon_ls = norm(S*w_ls - p_d,2);
%%
w_ls002 = pinv(S_002)*p_d_002;      % This solution gives the LS error
epsilon_ls002 = norm(S_002*w_ls002 - p_d_002,2);
 
w_ls005 = pinv(S_005)*p_d_005;      % This solution gives the LS error
epsilon_ls005 = norm(S_005*w_ls005 - p_d_005,2);

w_ls008 = pinv(S_008)*p_d_008;      % This solution gives the LS error
epsilon_ls008 = norm(S_008*w_ls008 - p_d_008,2);

%% LS solution with cvx minimization
% no relaxation
cvx_begin
    variable w_ls_cvx(N) complex
    minimize(norm(S*w_ls_cvx - p_d))
cvx_end
epsilon_ls_cvx = norm(S*w_ls_cvx - p_d,2);
% 0.02 relaxation
cvx_begin
    variable w_ls_cvx002(N) complex
    minimize(norm(S_002*w_ls_cvx002 - p_d_002))
cvx_end
epsilon_ls_cvx002 = norm(S_002*w_ls_cvx002 - p_d_002,2);
% 0.05 relaxation
cvx_begin
    variable w_ls_cvx005(N) complex
    minimize(norm(0.01*S_005*w_ls_cvx005 - p_d_005))
cvx_end
epsilon_ls_cvx005 = norm(S_005*w_ls_cvx005 - p_d_005,2);
% 0.08 relaxation
cvx_begin
    variable w_ls_cvx008(N) complex
    minimize(norm(S_008*w_ls_cvx008 - p_d_008))
cvx_end
epsilon_ls_cvx008 = norm(S_008*w_ls_cvx008 - p_d_008,2);
% The above minimization has a problem with scaling. So a factor 1/1000
% are used in the front of the S matrix. 

%% LS solution with cvx minimization with improved scaling
% no relaxation
cvx_begin
    variable w_ls_cvx(N) complex
    minimize(norm(S*w_ls_cvx/L - p_d))
cvx_end
epsilon_ls_cvx_L = cvx_optval;
%
cvx_begin
    variable w_ls_cvx_Lsquared(N) complex
    minimize(norm(S*w_ls_cvx_Lsquared/L^2 - p_d))
cvx_end
epsilon_ls_cvx_Lsquared = cvx_optval;

%
% 0.02 relaxation
cvx_begin
    variable w_ls_cvx002(N) complex
    minimize(norm(S_002*w_ls_cvx002/L^2 - p_d_002))
cvx_end
epsilon_ls_cvx002_L = cvx_optval;
%
cvx_begin
    variable w_ls_cvx002_Lsquared(N) complex
    minimize(norm(S_002*w_ls_cvx002_Lsquared/L^2 - p_d_002))
cvx_end
epsilon_ls_cvx002_Lsquared = cvx_optval;

%
% 0.05 relaxation
cvx_begin
    variable w_ls_cvx005(N) complex
    minimize(norm(0.01*S_005*w_ls_cvx005/L - p_d_005))
cvx_end
epsilon_ls_cvx005_L = cvx_optval;
%
cvx_begin
    variable w_ls_cvx005_Lsquared(N) complex
    minimize(norm(S_005*w_ls_cvx005_Lsquared/L^2 - p_d_005))
cvx_end
epsilon_ls_cvx005_Lsquared = cvx_optval;

%
% 0.08 relaxation
cvx_begin
    variable w_ls_cvx008(N) complex
    minimize(norm(S_008*w_ls_cvx008/L - p_d_008))
cvx_end
epsilon_ls_cvx008_L = cvx_optval;
%
cvx_begin
    variable w_ls_cvx008_Lsquared(N) complex
    minimize(norm(S_008*w_ls_cvx008_Lsquared/L^2 - p_d_008))
cvx_end
epsilon_ls_cvx008_Lsquared = cvx_optval;
%%
% Compare transition relaxation with 0 degree scaled with 1/L and 1/L^2
figure
subplot(2,2,1)
stem(abs(w_ls_cvx/L));
title('No relaxation scaled with 1/L')
subplot(2,2,2)
plot(angle,20*log10(abs(S*w_ls_cvx/L)))
title('No relaxation scaled with 1/L')

subplot(2,2,3)
stem(abs(w_ls_cvx_Lsquared/L^2));
title('No relaxation scaled with 1/L^2')
subplot(2,2,4)
plot(angle,20*log10(abs(S*w_ls_cvx_Lsquared/L^2)))
title('No relaxation scaled with 1/L^2')

% Compare transition relaxation with 002 degree scaled with 1/L and 1/L^2
figure
subplot(2,2,1)
stem(abs(w_ls_cvx002/L));
title('002 relaxation scaled with 1/L')
subplot(2,2,2)
plot(angle,20*log10(abs(S*w_ls_cvx002/L)))
title('002 relaxation scaled with 1/L')

subplot(2,2,3)
stem(abs(w_ls_cvx002_Lsquared/L^2));
title('002 relaxation scaled with 1/L^2')
subplot(2,2,4)
plot(angle,20*log10(abs(S*w_ls_cvx002_Lsquared/L^2)))
title('002 relaxation scaled with 1/L^2')

% Compare transition relaxation with 005 degree scaled with 1/L and 1/L^2
figure
subplot(2,2,1)
stem(abs(w_ls_cvx005/L));
title('005 relaxation scaled with 1/L')
subplot(2,2,2)
plot(angle,20*log10(abs(S*w_ls_cvx005/L)))
title('005 relaxation scaled with 1/L')

subplot(2,2,3)
stem(abs(w_ls_cvx005_Lsquared/L^2));
title('005 relaxation scaled with 1/L^2')
subplot(2,2,4)
plot(angle,20*log10(abs(S*w_ls_cvx005_Lsquared/L^2)))
title('005 relaxation scaled with 1/L^2')

% Compare transition relaxation with 008 degree scaled with 1/L and 1/L^2
figure
subplot(2,2,1)
stem(abs(w_ls_cvx008/L));
title('008 relaxation scaled with 1/L')
subplot(2,2,2)
plot(angle,20*log10(abs(S*w_ls_cvx008/L)))
title('008 relaxation scaled with 1/L')

subplot(2,2,3)
stem(abs(w_ls_cvx008_Lsquared/L^2));
title('008 relaxation scaled with 1/L^2')
subplot(2,2,4)
plot(angle,20*log10(abs(S*w_ls_cvx008_Lsquared/L^2)))
title('008 relaxation scaled with 1/L^2')

%% Exam the LS solution with cvx
figure
subplot(4,2,1)
stem(abs(w_ls_cvx/L));
subplot(4,2,2)
plot(angle,20*log10(abs(S*w_ls_cvx/L)))
subplot(4,2,3)
stem(abs(w_ls_cvx002/L));
subplot(4,2,4)
plot(angle,20*log10(abs(S*w_ls_cvx002/L)))
subplot(4,2,5)
stem(abs(w_ls_cvx005/L));
subplot(4,2,6)
plot(angle,20*log10(abs(S*w_ls_cvx005/L)))
subplot(4,2,7)
stem(abs(w_ls_cvx008/L));
subplot(4,2,8)
plot(angle,20*log10(abs(S*w_ls_cvx008/L)))

%% Exam the LS solution with pinv()
figure
subplot(4,2,1)
stem(abs(w_ls));
subplot(4,2,2)
plot(angle,20*log10(abs(S*w_ls)))
subplot(4,2,3)
stem(abs(w_ls002));
subplot(4,2,4)
plot(angle,20*log10(abs(S*w_ls002)))
subplot(4,2,5)
stem(abs(w_ls005));
subplot(4,2,6)
plot(angle,20*log10(abs(S*w_ls005)))
subplot(4,2,7)
stem(abs(w_ls008));
subplot(4,2,8)
plot(angle,20*log10(abs(S*w_ls008)))


% 0.0054

%% Solve the LASSO with different epsilon choice
epsilon_005_1 = linspace(epsilon_ls005,epsilon_max,25);
W_005_1 = zeros(N,size(epsilon_005_1,2));

%%
for i = 1:1:size(epsilon_005_1,2)
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S*w - p_d,2) <= epsilon_005_1(i);
cvx_end
W_005_1(:,i) = w;
end





