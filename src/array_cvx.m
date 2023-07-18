%% Array thinning cvx
close all
clc
clear
%% Parameters
c = physconst('LightSpeed');
lambda = 1550e-9;
fc = c/lambda;
k = 2*pi/lambda;
N = 300;                % to have a aperture of 1020 lambda
d = 0.5*lambda;
L = 3000;             % The resolution is at least 0.035 degree
leftend = -90;
rightend = 90;
Res = (rightend - leftend)/(L - 1);
angle = linspace(leftend,rightend,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';

% elementPos = (-(N-1)*d/2:d:(N-1)*d/2)';
% ang = 0;        % target angle
% w = ones(N,1);

%% S matrix
% S = transpose(exp(1i*k*x*sin(theta')));
S = exp(transpose(1i*k*x*sin(theta')));   % alternative
S2 = exp(1i*k*sin(theta)*x');
% S and S2 are suppose to be the same???
% But why it's not zero when subtrating them
%%
% AF_opt = AF(elementPos,-90,90,30000,k,ang,true,w);
%% Set the desired beam pattern
e = 1;                  % The error margin
% sll = -15;              % in dB
desired_pattern = upperboundgen(leftend,rightend,L);        % The desired pattern in 
%%
plot(angle,desired_pattern)

%%
cvx_begin
    variable w1(N) complex
    minimize(norm(w1,1))
    subject to
        % 20*log10(abs(S*w/nnz(w)))               % N is the # of active elements
        % 20*log(abs(S*w/nnz(w))/log(10))
        % not good, having problem here.......
        norm(S2*w1 - desired_pattern,2) <= 3;       % problem, not in dB
cvx_end

%%
plot(angle,20*log10(abs(S2*w)))

figure;
plot(abs(w1))
% choosing a value as threshold, counting the # of nonzero element
nnz(max(abs(w1)-0.05,0))


%%
cvx_begin
    variable w2(N) complex
    minimize(norm(w2,1))
    subject to
        % 20*log10(abs(S*w/nnz(w)))               % N is the # of active elements
        % 20*log(abs(S*w/nnz(w))/log(10))
        % not good, having problem here.......
        norm(S2*w2 - desired_pattern,2) <= 8;       % problem, not in dB
cvx_end
%%
plot(angle,20*log10(abs(S2*w2)))

figure;
plot(abs(w2))
% choosing a value as threshold, counting the # of nonzero element
nnz(max(abs(w2)-0.05,0))


%%
function reference_pattern = upperboundgen(left,right,L)   % only symmetry rect 
% sll in dB, 
% left and right indicate the range of the main lobe
% L is the sampling number with the whole FOV
    % sll_value = 10^(sll/20);
    
    syms y(x)
    y(x) = piecewise(abs(x) <= 30,1,0);  
    x = linspace(left,right,L)';
    reference_pattern = double(y(x));
    
end