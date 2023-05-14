close all
clc
clear
%% Try to find a better solution than the paper
% Use a small number of element to verify the priciple
c = physconst('LightSpeed');
lambda = 1550e-9;
fc = c/lambda;
N = 200;
d = 0.5*lambda;
Res = 3000;            % The resolution is at least 0.035 degree
leftend = -90;
rightend = 90;
k = 2*pi/lambda;
elementPos = (-(N-1)*d/2:d:(N-1)*d/2)';
ang = 0;        % target angle
% AFdB = AF(elementPos,leftend,rightend,Res,k,ang,true);
%% Genetic algorithm below
% GA's parameters
iter = 1000;                 % # of iterations
numG = 200;          % # of population
prvS = rng(2023);           % set random seed
% Step 3 Generating the random list
w0 = rand(N/2,numG)>0.5;     % initial generation with 200 candidates(symmetry)
% template candidate
% wtempr = w0(:,100);
% wtempl = flipud(wtempr);
% wtemp = [wtempl;wtempr];
% AF(elementPos,leftend,rightend,Res,k,ang,true,wtemp);

w0(:,99) = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 0 0 1 1 1 1 1 1 1 0 1 1 1 1 1 0 1 0 0 0 0 1 0 1 0 1 1 1 0 0 1 0 0 1 0 1 1 1 0 0 0 1 0 1 0 1 1 0 1];

%% GA
wm = w0;
% wl = flipud(wr);
% w = [wl;wr];                % introduce symmetry
AF_ = zeros(Res,numG);
pks = zeros(2,numG);
for m = 1:1:iter
    % 
    wl = flipud(wm);
    w = [wl;wm]; 
    % Calculate the AF of each eleemnts
    for i = 1:1:numG
        AF_(:,i) = AF(elementPos,leftend,rightend,Res,k,ang,false,w(:,i));
    end
    % AF_ has a row as a whole AF, # of rows is the numebr of populations
    for i  = 1:1:numG
        pks(:,i) = findpeaks(AF_(:,i),'NPeaks',2,'SortStr','descend');
    end
    sll = 0 - pks(2,:);
    % Sort the resulting sidelobe level
    [~,idx] = sort(sll,2,'descend');
    % (Selection) Discard half of the generation that gets the lower score 
    wm = wm(:,[idx(1:numG/2) idx(1:numG/2)]);
    % Paring and mating
    
    % Mutate 
    mutIdx = randi(N/2,1,3);
    wm(mutIdx,numG/2+1:numG) = not(wm(mutIdx,numG/2+1:numG));
end

%%
woptr = wm(:,1);
woptl = flipud(woptr);
wopt = [woptl;woptr];      
AF(elementPos,-90,90,Res,k,ang,true,wopt);
% figure
% plot(elementPos_opt,zeros(size(elementPos_opt)),'o')
% AF(elementPos,-90,90,Res,k,N,ang,true);

%%
disp_fr(wm,N,numG)