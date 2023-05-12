close all
clc
clear
%% Try to reshow the result in the literature
% Use a small number of element to verify the priciple
c = physconst('LightSpeed');
lambda = 1550e-9;
fc = c/lambda;
% N = 200;
N = 200;
d = 0.5*lambda;
Res = 3000;            % The resolution is at least 0.035 degree
leftend = -90;
rightend = 90;
k = 2*pi/lambda;
elementPos = (-(N-1)*d/2:d:(N-1)*d/2)';
ang = 0;        % target angle

% w1 = ones(size(elementPos));
AFdB = AF(elementPos,leftend,rightend,Res,k,ang,true);
%% Genetic algorithm below
prvS = rng(2023);           % set random seed
% Step 3 Generating the random list
w0 = rand(N,1,200)>0.5;     % initial generation with 200 candidates
% template candidate
wtemp = w0(:,:,100);
% wtemp(1:100) = 1;
% elementPostemp = wtemp.*elementPos;
AF(elementPos,leftend,rightend,Res,k,ang,true,wtemp);
% GA's parameters
iter = 100;                  % # of iterations
numG = size(w0,3);       % # of population
w = w0;
% elementPos_ = w.*elementPos;
AF_ = zeros(Res,numG);
pks = zeros(2,numG);
for m = 1:1:iter
    % Calculate the AF of each eleemnts
    for i = 1:1:numG
        AF_(:,i) = AF(elementPos,leftend,rightend,Res,k,ang,false,w(:,:,i));
    end
    % AF_ has a row as a whole AF, # of rows is the numebr of populations
    for i  = 1:1:numG
        pks(:,i) = findpeaks(AF_(:,i),'NPeaks',2,'SortStr','descend');
    end
    sll = 0 - pks(2,:);
    % Sort the resulting sidelobe level
    [~,idx] = sort(sll,2,'descend');
    % (Selection) Discard half of the generation that gets the lower score 
    w = w(:,:,[idx(1:numG/2) idx(1:numG/2)]);
    % Paring and mating

    % Mutate 
    mutIdx = randi(N,1,3);
    w(mutIdx,:,numG/2+1:numG) = not(w(mutIdx,:,numG/2+1:numG));
end
%%
wopt = w(:,:,1);
AF(elementPos,-90,90,Res,k,ang,true,wopt);
% figure
% plot(elementPos_opt,zeros(size(elementPos_opt)),'o')

%% Continue GA
% GA's parameters
iter = 50;                  % # of iterations
for m = 1:1:iter
    % seems lack of a elementPos_ = w.*elementPos; updating
    % Calculate the AF of each eleemnts
    for i = 1:1:numG
        AF_(:,i) = AF(elementPos,leftend,rightend,Res,k,ang,false,w(:,:,i));
    end
    % AF_ has a row as a whole AF, # of rows is the numebr of populations
    for i  = 1:1:numG
        pks(:,i) = findpeaks(AF_(:,i),'NPeaks',2,'SortStr','descend');
    end
    sll = 0 - pks(2,:);
    % Sort the resulting sidelobe level
    [~,idx] = sort(sll,2,'descend');
    % (Selection) Discard half of the generation that gets the lower score 
    w = w(:,:,[idx(1:numG/2) idx(1:numG/2)]);
    % Paring and mating

    % Mutate 
    mutIdx = randi(N,1,3);
    w(mutIdx,:,numG/2+1:numG) = not(w(mutIdx,:,numG/2+1:numG));
end

%%
wopt = w(:,:,1);
AF(elementPos,-90,90,Res,k,ang,true,wopt);
% figure
% plot(elementPos_opt,zeros(size(elementPos_opt)),'o')
% AF(elementPos,-90,90,Res,k,N,ang,true);