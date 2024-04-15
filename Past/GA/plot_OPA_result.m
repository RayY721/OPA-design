close all
clc
clear
%% Plot the result for OPA
c = physconst('LightSpeed');
lambda = 1550e-9;
fc = c/lambda;
N = 600;                % to have a aperture of 1020 lambda
d = 1.7*lambda;         
Res = 3000;             % The resolution is at least 0.035 degree
leftend = -20;
rightend = 20;
k = 2*pi/lambda;
elementPos = (-(N-1)*d/2:d:(N-1)*d/2)';
ang = 0;        % target angle

load GA_OPA1.mat;

%% 
woptr = wr(:,1);
woptl = flipud(woptr);
wopt = [woptl;woptr];      
AF_opt = AF(elementPos,-90,90,30000,k,ang,true,wopt);
title('The AF over -90 and 90 degree')

AF_opt_narrow = AF(elementPos,leftend,rightend,Res,k,ang,true,wopt);
title('The AF over the region of interest (-20,20)')
figure
findpeaks(AF_opt_narrow,'NPeaks',2,'SortStr','descend')