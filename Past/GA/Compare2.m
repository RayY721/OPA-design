% Compare the result from 'GA_intensively_mutation' and the paper. 
lambda = 1550e-9;
d = 0.5*lambda;
Res = 3000;
k = 2*pi/lambda;
N = 200;
elementPos = (-(N-1)*d/2:d:(N-1)*d/2)';
ang = 0;        % target angle

load wopt.mat;

%%
wr = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 0 0 1 1 1 1 1 1 1 0 1 1 1 1 1 0 1 0 0 0 0 1 0 1 0 1 1 1 0 0 1 0 0 1 0 1 1 1 0 0 0 1 0 1 0 1 1 0 1]';

wl = flipud(wr);
w_reprod = [wl;wr];
AF_reprod = AF(elementPos,-90,90,Res,k,ang,true,w_reprod);
ylim([-40 0])
title('The result from reference')
figure
findpeaks(AF_reprod,'NPeaks',2,'SortStr','descend')
title('The max side lobe of reference result')

woptl = flipud(woptr);
wopt = [woptl;woptr];      
AF_opt = AF(elementPos,-90,90,Res,k,ang,true,wopt);
ylim([-40 0])
title('The result from mutation only search')
figure
findpeaks(AF_opt,'NPeaks',2,'SortStr','descend')
title('The max side lobe of result from mutation only search')