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
L = 700;             % The resolution is at least 0.035 degree
win.leftend = -18;
win.rightend = 18;
Res = (win.rightend - win.leftend)/(L - 1);
angle = linspace(win.leftend,win.rightend,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
%% Set the desired beam pattern
left_win = -0.04;
right_win = 0.04;
desired_pattern = upperboundgen(win.leftend,win.rightend,[left_win right_win],L);       % The desired pattern
figure
plot(angle,desired_pattern)
%% Modify the S matrix and reference pattern    % have no effect when L = 700/701
loose_range = 0.005;    % in degree
S_new = adjustSmatrix(S,[win.leftend win.rightend],[left_win right_win],L,loose_range);
desired_pattern_modified = adjustRefpattern(desired_pattern,[win.leftend win.rightend],[left_win right_win],L,loose_range);

%%
cvx_begin
    variable w(N) complex
    minimize(norm(w,1))
    subject to
        norm(S_new*w - desired_pattern_modified,2) <= 0.21;       % problem, not in dB
        
cvx_end
%% evaluating result with a higher resolution
L = 3301;
win.leftend = -18;
win.rightend = 18;
Res = (win.rightend - win.leftend)/(L - 1);
angle = linspace(win.leftend,win.rightend,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
figure
subplot(2,1,1)
plot(angle,20*log10(abs(S*w)))
title("beampattern of w with L =" + num2str(L))
subplot(2,1,2)
stem(abs(w))
%%
L = 700;             % The resolution is at least 0.035 degree
win.leftend = -18;
win.rightend = 18;
Res = (win.rightend - win.leftend)/(L - 1);
angle = linspace(win.leftend,win.rightend,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
desired_pattern = upperboundgen(win.leftend,win.rightend,[left_win right_win],L);   
figure
stem(abs(S*w - desired_pattern))
title('The difference between generated pattern and reference pattern when L=700')
%%
L = 701;
win.leftend = -18;
win.rightend = 18;
Res = (win.rightend - win.leftend)/(L - 1);
angle = linspace(win.leftend,win.rightend,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
desired_pattern = upperboundgen(win.leftend,win.rightend,[left_win right_win],L); 
figure
stem(abs(S*w - desired_pattern))
title('The difference between generated pattern and reference pattern when L=701')



%% Parameters
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;
N = 1000;                % to have a aperture of 1020 lambda
d = param.lambda;
L = 700;             % The resolution is at least 0.035 degree
win.leftend = -18;
win.rightend = 18;
Res = (win.rightend - win.leftend)/(L - 1);
angle = linspace(win.leftend,win.rightend,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
% Set the desired beam pattern
left_win = -0.04;
right_win = 0.04;
desired_pattern = upperboundgen(win.leftend,win.rightend,[left_win right_win],L);       % The desired pattern
figure
plot(angle,desired_pattern)
% Modify the S matrix and reference pattern    % have no effect when L = 700/701
loose_range = 0.005;    % in degree
S_new = adjustSmatrix(S,[win.leftend win.rightend],[left_win right_win],L,loose_range);
desired_pattern_modified = adjustRefpattern(desired_pattern,[win.leftend win.rightend],[left_win right_win],L,loose_range);
%
cvx_begin
    variable w2(N) complex
    minimize(norm(w2,1))
    subject to
        norm(S_new*w2 - desired_pattern_modified,2) <= 0.61;       % problem, not in dB
        
cvx_end

%%
L = 3301;
win.leftend = -18;
win.rightend = 18;
Res = (win.rightend - win.leftend)/(L - 1);
angle = linspace(win.leftend,win.rightend,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
figure
subplot(2,1,1)
plot(angle,20*log10(abs(S*w2)))
title("beampattern of w2 with L =" + num2str(L))
subplot(2,1,2)
stem(abs(w2))
% desired_pattern = upperboundgen(win.leftend,win.rightend,[left_win right_win],L); 

%%
figure
subplot(2,1,1)
semilogy(flip(sort(abs(w))))
subplot(2,1,2)
semilogy(flip(sort(abs(w2))))