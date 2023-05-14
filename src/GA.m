clc
clear
close all
%% The Genetic Algorithm

c = physconst('LightSpeed');
lambda = 1550e-9;
fc = c/lambda;
% for minimum spacing of 2500nm, which is about 1.7 lambda. 
% The aperture is about 1020 lambda. 
% So the fully positioned array, N = 600
N = 600;
d = 1.7*lambda;

Res = 3000;            % The resolution is at least 0.035 degree
leftend = -20;
rightend = 20;
k = 2*pi/lambda;
elementPos = (-(N-1)*d/2:d:(N-1)*d/2)';

% clf;
% wplot = helperThinnedArrayComparison(elementPos,fc,c);
ang = 0;        % target angle
AFdB = AF(elementPos,leftend,rightend,Res,k,N,ang,true);

%% genetic algorithm below
prvS = rng(2023);           % set random seed
% Step 3 Generating the random list
w0 = rand(N,1,200)>0.5;     % initial generation with 200 candidates
% template candidate
wtemp = w0(:,:,100);
% wtemp(1:100) = 1;
elementPostemp = wtemp.*elementPos;
AF(elementPostemp,leftend,rightend,Res,k,N,ang,true);
% GA's parameters
iter = 1000;                  % # of iterations
numGene = size(w0,3);       % # of population
w = w0;
elementPos_ = w0.*elementPos;
AF_ = zeros(Res,numGene);
pks = zeros(2,numGene);
for m = 1:1:iter
    % Calculate the AF of each eleemnts
    for i = 1:1:numGene
        AF_(:,i) = AF(elementPos_(:,:,i),leftend,rightend,Res,k,N,ang,false);
    end
    % AF_ has a row as a whole AF, # of rows is the numebr of populations
    for i  = 1:1:numGene
        pks(:,i) = findpeaks(AF_(:,i),'NPeaks',2,'SortStr','descend');
    end
    sll = 0 - pks(2,:);
    % Sort the resulting sidelobe level
    [~,idx] = sort(sll,2,'descend');
    % Discard half of the generation that gets the lower score 
    w = w(:,:,[idx(1:numGene/2) idx(1:numGene/2)]);
    % Paring and mating

    % Mutate 
    mutIdx = randi(N,1,3);
    w(mutIdx,:,numGene/2+1:numGene) = not(w(mutIdx,:,numGene/2+1:numGene));

    % updating the elements' position
    elementPos_ = w.*elementPos;
end
%%
wopt = w(:,:,1);
elementPos_opt = wopt.*elementPos;
AF(elementPos_opt,-90,90,Res,k,N,ang,true);
% figure
% plot(elementPos_opt,zeros(size(elementPos_opt)),'o')
% AF(elementPos,-90,90,Res,k,N,ang,true);

%% Continue GA
% GA's parameters
iter = 500;                  % # of iterations
for m = 1:1:iter
    % seems lack of a elementPos_ = w.*elementPos; updating
    % Calculate the AF of each eleemnts
    for i = 1:1:numGene
        AF_(:,i) = AF(elementPos_(:,:,i),leftend,rightend,Res,k,N,ang,false);
    end
    % AF_ has a row as a whole AF, # of rows is the numebr of populations
    for i  = 1:1:numGene
        pks(:,i) = findpeaks(AF_(:,i),'NPeaks',2,'SortStr','descend');
    end
    sll = 0 - pks(2,:);
    % Sort the resulting sidelobe level
    [~,idx] = sort(sll,2,'descend');
    % (Selection) Discard half of the generation that gets the lower score 
    w = w(:,:,[idx(1:numGene/2) idx(1:numGene/2)]);
    % Paring and mating

    % Mutate 
    mutIdx = randi(N,1,3);
    w(mutIdx,:,numGene/2+1:numGene) = not(w(mutIdx,:,numGene/2+1:numGene));
    % updating the elements' position
    elementPos_ = w.*elementPos;
end
%%
wopt = w(:,:,1);
AF(elementPos_opt,-90,90,Res,k,N,ang,true);
%%
AF(elementPostemp,-90,90,Res,k,N,ang,true);