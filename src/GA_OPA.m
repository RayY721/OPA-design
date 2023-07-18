%% GA for OPA design
close all
clc
clear
%%
c = physconst('LightSpeed');
lambda = 1550e-9;
fc = c/lambda;
N = 50;                % to have a aperture of 1020 lambda
d = 1.7*lambda;         
Res = 3000;             % The resolution is at least 0.035 degree
leftend = -90;
rightend = 90;
k = 2*pi/lambda;
elementPos = (-(N-1)*d/2:d:(N-1)*d/2)';
ang = 0;        % target angle

%% GA's parameters
iter = 1000;
numG = 400;
prvS = rng(2023);           % set random seed
w = rand(N/2,numG)>0.5;
wr = GA(numG,Res,iter,leftend,rightend,k,ang,elementPos,w);

%%
disp_fr(wr,N,numG)
%%
woptr = wr(:,1);
woptl = flipud(woptr);
wopt = [woptl;woptr];      
AF_opt = AF(elementPos,-90,90,30000,k,ang,true,wopt);
% findpeaks(AF_opt,'NPeaks',2,'SortStr','descend')
%%
w = ones(N,1);
AF = AF(elementPos,-90,90,30000,k,ang,true,w);