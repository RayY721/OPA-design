%% OMP for joint grid method
close all
clc
clear
% Parameters
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;
N = 10;                % to have a aperture of 1020 lambda
d = 0.5*param.lambda;
L = 500;             % now Lz = 700 or 900 doesn't work
win.leftend = -90;
win.rightend = 90;
Res = (win.rightend - win.leftend)/(L - 1);
angle = linspace(win.leftend,win.rightend,L)';  % in degree
theta = angle*pi/180;
u = linspace(sin(win.leftend*pi/180),sin(win.rightend*pi/180),L)';
x = (-(N-1)*d/2:d:(N-1)*d/2)';
% x = (0:d:(N-1)*d)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
S_u = exp(1i*param.k*u*x');
% Discretize the unit circle
M = 100;   
psi = linspace(0,2*pi,M+1);   % actual input should be desired number plus 1
psi = psi(1:end-1)';
b = exp(1i*psi);
% polarplot(angle(b),abs(b),"o")
% Generate the model y = alpha*Ax
A = kron(S,b.');
% Set the desired beam pattern (Loose than our requirement) 
left_win = -2;
right_win = 2;
% left_win = -0.04;
% right_win = 0.04;

desired_pattern = upperboundgen(win.leftend,win.rightend,[left_win right_win],L);       % The desired pattern
figure
plot(angle,desired_pattern)
xlabel('Angle (degrees)');
ylabel('Amplitude')
title('Ideal beam pattern')
%% inspect the "sensing matrix"
[X,Y] = meshgrid(1:1000,1:500);
surf(X,Y,real(A));
%% Modify the S matrix and reference pattern
loose_range = 0.02;    % in degree   (0.014 is not feasible)
% loose_range = 0.2;

desired_pattern_modified = adjustrefpattern(desired_pattern,[win.leftend win.rightend],[left_win right_win],L,loose_range);

S_new = adjustSmatrix(S,[win.leftend win.rightend],[left_win right_win],L,loose_range);

%% Coherence between column of A and the y
% coherence_vector = zeros(size(A,2),1);
% for i = size(A,2)
%     coherence_vector(i)=abs(A(:,i)'*desired_pattern);
% end
support = [];
deleted = [];
r = desired_pattern;
z = zeros(size(A,2),1);
%%
% for i = 1:1:size(A,2)

coherence_vector=abs(A'*r);
coherence_vector(deleted) = 0;
[val,idx] = max(coherence_vector);
support = [support, idx];
% very unlikely need the following two lines, as what matters is the
% amplitude of the array factor (AF). The difference should also be in the
% amplitude, rather than the true complex domain. 
% b_ = [real(b);imag(b)];
% a_ = [real(desired_pattern;imag(desired_pattern)];

% alpha = pinv(abs(A(:,support)))*desired_pattern;    % what fuck? why this is not a scaler???????
z(support) = 1;
% r = desired_pattern - alpha(1)*abs(A(:,support)*z(support));
r = desired_pattern - abs(A(:,support)*z(support));
deleted = [deleted,(ceil(idx/M)-1)*M+1:ceil(idx/M)*M]; % ceil(idx/M) is the No. of groups, nth group 
% Need to delete No.ceil(idx/M) to No.floor(idx/M)+M


subplot(2,1,1)
plot(angle,abs(r))
subplot(2,1,2)
% plot(angle,20*log10(alpha(1)*abs(A(:,support)*z(support))))
plot(angle,20*log10(abs(A(:,support)*z(support))))
% end
%%
figure
subplot(5,1,1)
plot(angle,real(A(:,25)))
subplot(5,1,2)
plot(angle,real(A(:,24)))
subplot(5,1,3)
plot(angle,real(A(:,23)))
subplot(5,1,4)
plot(angle,real(A(:,22)))
subplot(5,1,5)
plot(angle,abs(A(:,21)))
%% Function OMP
% Given "measurement" y, "CS matrix" A, "sparsity" k
% function [z,err] = OMP(y,A,epsilon)
function z = OMP(y,A,epsilon)
N = max(size(A));
support=[];
r=y;  
z = zeros(N,1);
% prod = zeros(N,1);
for times = 1:1:10000
    prod = abs(A'*r);
    [val,pos] = max(prod);
    support = [support, pos];
    % LS
    z_lambda = pinv(A(:,support))*y;
    r = y - A(:,support)*z_lambda;
%     err_(times) = norm(r,2);
    err_ = norm(r,2)^2;
    if err_<epsilon*10 %
        err_
            break; 
    end
end
% err = err_';
z(support) = z_lambda;
end

function [x_estimate, support] = omp(A, y, k)
    % Orthogonal Matching Pursuit algorithm
    % Inputs:
    %   A: Measurement matrix
    %   y: Measurement vector
    %   k: Desired sparsity (number of non-zero elements in x)
    % Outputs:
    %   x_estimate: Estimated sparse signal
    %   support: Indices of non-zero elements in x

    [N, M] = size(A);
    x_estimate = zeros(M, 1);
    r = y; % Residual

    support = []; % Indices of non-zero elements in x

    for iter = 1:k
        % Compute inner products and find the index corresponding to the maximum inner product
        inner_products = abs(A' * r);
        [~, idx] = max(inner_products);

        % Update support set
        support = [support; idx];

        % Solve least squares problem using the current support set
        x_temp = A(:, support)\y;

        % Update the estimate of x
        x_estimate(support) = x_temp;

        % Update the residual
        r = y - A * x_estimate;

        % Check for convergence (small residual)
        if norm(r) < 1e-6
            break;
        end
    end
end

function [x_estimate, support] = omp_binary_block_structure(A, y, N, M)
    % OMP algorithm with binary sparse vector and block structure
    % Inputs:
    %   A: Measurement matrix
    %   y: Measurement vector
    %   N: Number of blocks
    %   M: Number of entries in each block
    % Outputs:
    %   x_estimate: Estimated sparse binary vector with block structure
    %   support: Indices of selected entries in x_estimate

    [P, Q] = size(A);
    x_estimate = zeros(Q, 1);
    r = y; % Residual

    support = []; % Indices of selected entries in x_estimate

    for block = 1:N
        % Indices corresponding to the current block
        block_indices = (block - 1) * M + 1 : block * M;

        % Compute inner products and find the index corresponding to the maximum inner product within the current block
        inner_products = abs(A' * r);
        inner_products(block_indices) = 0; % Exclude already selected indices
        [~, idx_within_block] = max(inner_products(block_indices));

        % Update support set
        support = [support; block_indices(idx_within_block)];

        % Solve least squares problem using the current support set
        x_temp = A(:, support)\y;

        % Update the estimate of x
        x_estimate(support) = x_temp;

        % Update the residual
        r = y - A * x_estimate;

        % Check for convergence (small residual)
        if norm(r) < 1e-6
            break;
        end
    end
end
