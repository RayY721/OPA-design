%% Compare the Array pattern in Haupt's paper and the cvx I formulated
close all
clc
clear
%%
c = physconst('LightSpeed');
lambda = 1550e-9;
fc = c/lambda;
k = 2*pi/lambda;
N = 200;
d = 0.5*lambda;
L = 3000;            % The resolution is at least 0.035 degree
leftend = -90;
rightend = 90;
Res = (rightend - leftend)/(L - 1);
% angle = linspace(leftend,rightend,L)';  % in degree
% theta = angle*pi/180;
angle = -90:0.01:90;
theta = angle'*pi/180;
% x = (0:d:(N-1)*d)';
x = (-(N-1)*d/2:d:(N-1)*d/2)';
WAVE.TAR_THETA = 0;      % target direction
WAVE.HALF_BEAMWIDTH =1 ;       % expected main lobe width

%% GA' result
wr = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 0 0 1 1 1 1 1 1 1 0 1 1 1 1 1 0 1 0 0 0 0 1 0 1 0 1 1 1 0 0 1 0 0 1 0 1 1 1 0 0 0 1 0 1 0 1 1 0 1]';
wl = flipud(wr);
w = [wl;wr]; 
%%
S = exp(transpose(1i*k*x*sin(theta')));   % alternative
S2 = exp(1i*k*sin(theta)*x');
%% 
plot(angle',20*log10(abs(S*w/nnz(w))))
ylim([-40 0])
figure
plot(angle',20*log10(abs(S2*w/nnz(w))))
ylim([-40 0])
%% 
Star = S(find(angle==WAVE.TAR_THETA),:);
ind = find(leftend<=angle&angle <= (WAVE.TAR_THETA-WAVE.HALF_BEAMWIDTH) | ...
           angle >= (WAVE.TAR_THETA+WAVE.HALF_BEAMWIDTH)&angle<=rightend );
Ss = S(ind,:);

%%
cvx_begin
    variable w(N) complex
    variable t
    minimize( norm(w,1) + t )
    subject to 
        Star*w == 1;
        abs(Ss*w) <= t
cvx_end
%%
plot(sin(theta),20*log10(abs(S*w)))