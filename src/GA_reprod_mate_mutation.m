close all
clc
clear
%%
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
%% Genetic algorithm below
% GA's parameters
iter = 3000;                 % # of iterations
numG = 400;          % # of population
prvS = rng(2023);           % set random seed
wrefr = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 0 0 1 1 1 1 1 1 1 0 1 1 1 1 1 0 1 0 0 0 0 1 0 1 0 1 1 1 0 0 1 0 0 1 0 1 1 1 0 0 0 1 0 1 0 1 1 0 1]';
% Step 3 Generating the random list
w0 = rand(N/2,numG)>0.5;     % initial generation with 200 candidates(symmetry)
w0(:,99) = wrefr;
%% GA
% inputs: wr, numG, Res, iter, leftend, rightend, k, ang, 
% output: wr
% some hyperparameters:
% way of crossover: single point crossover, mulipoint crossover and 
wr = w0;               
AF_ = zeros(Res,numG);
pks = zeros(2,numG);
for m = 1:1:iter
    % 
    wl = flipud(wr);            % introduce symmetry
    w = [wl;wr]; 
    % Calculate the AF of each elements
    for i = 1:1:numG
        AF_(:,i) = AF(elementPos,leftend,rightend,Res,k,ang,false,w(:,i));
    end
    % AF_ has a row as a whole AF, # of rows is the numebr of populations
    for i  = 1:1:numG
        pks(:,i) = findpeaks(AF_(:,i),'NPeaks',2,'SortStr','descend');
    end
    sll = 0 - pks(2,:);
    % (Sort) the resulting sidelobe level
    [~,idx] = sort(sll,2,'descend');
    % (Selection) Discard half of the generation that gets the lower score 
    wr_ = wr(:,idx(1:numG/2));
    % (Paring and mating)
    crossover_point = randi(N/2-1,numG/4,1);    % half # of numG to be paired, quarter of numG crossover point
    for n = 1:1:numG/4
        offsprints = [wr_(1:crossover_point(n),2*n-1),wr_(1:crossover_point(n),2*n);
            wr_(crossover_point(n)+1:end,2*n),wr_(crossover_point(n)+1:end,2*n-1)];
        wr_ = [wr_,offsprints];
    end
    
    % Mutate 
    mutIdx = randi(N/2,1,1);
    % mutG = randi(numG/2,1,1)+100;
    % wr_(mutIdx,mutG) = not(wr_(mutIdx,mutG));
    wr_(mutIdx,numG/2+1:numG) = not(wr_(mutIdx,numG/2+1:numG));
    wr = wr_;
end
%%
woptr = wr(:,1);
woptl = flipud(woptr);
wopt = [woptl;woptr];      
AF_opt = AF(elementPos,-90,90,Res,k,ang,true,wopt);
findpeaks(AF_opt,'NPeaks',2,'SortStr','descend')
ylim([-40 0])
%%
disp_fr(wr,N,numG)
%% Check if equals
vec_eq = [];
for i = 1:1:199     % change the iteration number to compare the element
    vec_eq = [vec_eq isequal(wr(:,i),wr(:,i))];
end
nnz(vec_eq)
%% 
wrefl = flipud(wrefr);
wref = [wrefl;wrefr];
AF_ref = AF(elementPos,-90,90,Res,k,ang,true,wref);
ylim([-40 0])
findpeaks(AF_ref,'NPeaks',2,'SortStr','descend')