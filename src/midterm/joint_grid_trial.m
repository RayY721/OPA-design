% 
close all
clc
clear

N = 5;
M = 5;
d = 0.5;    % in wavelength

L = 500;

angle = linspace(-90,90,L)';  % in degree
theta = angle*pi/180;

x = (-(N-1)*d/2:d:(N-1)*d/2)';

S = exp(1i*2*pi*sin(theta)*x');  

psi = linspace(0,2*pi,M+1);   % actual input should be desired number plus 1
psi = psi(1:end-1)';
b = exp(1i*psi);

A = kron(S,b.');

%%
[X,Y] = meshgrid(1:25,1:500);
surf(X,Y,real(A))

%%
w = ones(N,1);
w(1) = 0;
figure
plot(angle,20*log10(abs(S*ones(N,1))))
figure
plot(angle,abs(S*ones(N,1)))

%%
d_dot = abs(S*ones(N,1));
d_dot(d_dot <= 2.2) = 0;
figure
plot(angle,d_dot);
%%
[X,Y] = meshgrid(1:N^2,1:500);
surf(X,Y,real(A));
%%
support = [];
deleted = [];
r = d_dot;
z = zeros(size(A,2),1);
%%
coherence_vector=abs(A'*r);
coherence_vector(deleted) = 0;
[val,idx] = max(coherence_vector);
%%
support = [support, idx];
% very unlikely need the following two lines, as what matters is the
% amplitude of the array factor (AF). The difference should also be in the
% amplitude, rather than the true complex domain. 
% b_ = [real(b);imag(b)];
% a_ = [real(desired_pattern;imag(desired_pattern)];

% alpha = pinv(abs(A(:,support)))*desired_pattern;    % what fuck? why this is not a scaler???????
z(support) = 1;
% r = desired_pattern - alpha(1)*abs(A(:,support)*z(support));
r = d_dot - abs(A(:,support)*z(support));
deleted = [deleted,(ceil(idx/M)-1)*M+1:ceil(idx/M)*M]; % ceil(idx/M) is the No. of groups, nth group 
% Need to delete No.ceil(idx/M) to No.floor(idx/M)+M
figure
subplot(2,1,1)
plot(angle,abs(r))
subplot(2,1,2)
% plot(angle,20*log10(alpha(1)*abs(A(:,support)*z(support))))
% plot(angle,20*log10(abs(A(:,support)*z(support))))
plot(angle,abs(A(:,support)*z(support)))
% end
%% plot all possible combination
figure
% for i = 1:1:25
%     subplot(5,5,i)
%     plot(angle,)
% 
% end
