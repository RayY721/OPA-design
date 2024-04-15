%% The GA converged to a solution with filling rate of 76%
% This script compare the array factor found by the
% script "GA_reprod_mate_mutation.m" and the reference paper. 
close all
clc
clear

%% 
lambda = 1550e-9;
d = 0.5*lambda;
Res = 3000;
k = 2*pi/lambda;
N = 200;
elementPos = (-(N-1)*d/2:d:(N-1)*d/2)';
ang = 0;        % target angle

load wopt_converged.mat;

%%
wrefl = flipud(wrefr);
wref = [wrefl;wrefr];
AF_ref = AF(elementPos,-90,90,Res,k,ang,true,wref);
ylim([-40 0])
title('The result from reference')
figure
findpeaks(AF_ref,'NPeaks',2,'SortStr','descend')
title('The max side lobe of reference result')

woptl = flipud(woptr);
wopt = [woptl;woptr];      
AF_opt = AF(elementPos,-90,90,Res,k,ang,true,wopt);
ylim([-40 0])
title('The result from crossover and mutation')
figure
findpeaks(AF_opt,'NPeaks',2,'SortStr','descend')
title('The max side lobe of result from algo with crossover and mutation')