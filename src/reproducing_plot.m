% This script reproduces the result of the "Thinned Arrays Using Genetic
% Algorithms" published on IEEE Transactions on Antennas and Propagation,
% by Randy L. Haupt. 
% The code below plot the Fig.3(a), Fig.5(a) and Fig.5(b)
% By plotting the patterns, the function AF has been verfied to be correct.
close all
clc
clear
%%
c = physconst('LightSpeed');
lambda = 1550;
% fc = c/lambda;
% N = 200;
N = 200;
d = 0.5*lambda;
Res = 5000;     
leftend = -90;
rightend = 90;
k = 2*pi/lambda;
elementPos = (-(N-1)*d/2:d:(N-1)*d/2)';
ang = 0;        % target angle
%%
for i = 1:1:3
    switch i
        case 1
            w = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 0 0 1 1 1 1 1 1 1 0 1 1 1 1 1 0 1 0 0 0 0 1 0 1 0 1 1 1 0 0 1 0 0 1 0 1 1 1 0 0 0 1 0 1 0 1 1 0 1];
            ang = 0;  
        case 2
            w = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 0 1 1 1 1 1 0 0 1 1 1 1 1 0 0 0 1 0 0 1 1 1 1 1 0 0 0 0 1 0 0 1 1 0 0 0 1 1 0 1 0 1 1 0 1 1 0 1];
            ang = 0;  
        case 3
            w = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 0 1 1 1 1 1 0 0 1 1 1 1 1 0 0 0 1 0 0 1 1 1 1 1 0 0 0 0 1 0 0 1 1 0 0 0 1 1 0 1 0 1 1 0 1 1 0 1];
            ang = -30;
        otherwise
            disp('error occurs')
    end
    wr = w';
    wl = flipud(wr);
    w_reprod = [wl;wr];
    AF(elementPos,-90,90,Res,k,ang,true,w_reprod);
    ylim([-40 0])
    switch i
        case 1
            title('Fig.3(a)');
        case 2
            title('Fig.5(a)');
        case 3
            title('Fig.5(b)');
        otherwise
            disp('error occurs')
    end
    findpeaks(AF(elementPos,-90,90,Res,k,ang,true,w_reprod),'NPeaks',2,'SortStr','descend')
    switch i
        case 1
            title('Fig.3(a)');
        case 2
            title('Fig.5(a)');
        case 3
            title('Fig.5(b)');
        otherwise
            disp('error occurs')
    end
end
% profile on

% profreport