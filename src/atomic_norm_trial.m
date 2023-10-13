%% for fun, try out the atomic norm

cvx_begin sdp
    variable u(L) complex;
    variable t complex;
    variable y(L) complex;
    T = toeplitz(u);

    minimize real((trace(Tu)/L + trace(t)))
    subject to 
        norm(desired_pattern - y) <= 0.1;
        [t y';y T] >= 0;

cvx_end

% https://github.com/JAMES-YI/00_Separation-free-Super-resolution/blob/main/SANM_cvx.m


% function [rx,err]=SANM_cvx(ox,K,N,kk,stepsize,epsilon)
% %% Standard Atomic Norm Minimization
% 
% % solve the dual problem
% cvx_quiet true
% cvx_precision default
% 
% cvx_begin sdp
%   variable rx(N) complex;
%   variable u(N) complex;
%   variable s complex;
%   Q = toeplitz(u);
% 
%   minimize real(trace(Q)/N + trace(s))
%   subject to 
%         [s rx'; rx Q] >= 0,
%         rx(K)==ox(K);
% cvx_end
% 
% err = norm(rx-ox)/norm(ox);
% end