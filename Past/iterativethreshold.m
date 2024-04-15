%% Iterative hard thresholding pursuit function
function [cost,threshold,w_result,W,MSE,MSE_new,MSE_app] = iterativethreshold(w,S,lambda)

% [~,index_max] = max(abs(S*w));

N = length(w);
cost = zeros(N,1);
MSE = zeros(N,1);
MSE_new = zeros(N,1);
MSE_app = zeros(N,1);
% MSE_difference = zeros(N,1);
w_t = w;
[~, sorted_indices] = sort(abs(w),'ascend');
L = size(S,1);

W = zeros(N);

% S_normalized = S/abs(S(index_max,:)*w);

% S_normalized = S./N;

S_app = diag(ones(N,1));

for i = 1:1:N
    
    index = sorted_indices(i);
    w_t(index) = 0;
    W(:,i) = w_t;
    % S_normalized_t = S/abs(S(index_max,:)*w_t);
    % cost function evaluation
    % cost(i) = lambda*nnz(w_t)/N + (w - w_t)'*S'*S*(w - w_t)/L;
    cost(i) = lambda*nnz(w_t)/N + (S*w - S*w_t)'*(S*w - S*w_t)/L;
    % MSE(i) = (S*w - S*w_t)'*(S*w - S*w_t)/L;
    
    e = w - w_t;
    MSE_new(i) = (S*e)'*(S*e)/L;
    MSE_app(i) = e'*S_app*e/N;
    MSE(i) = (S*w - S*w_t)'*(S*w - S*w_t)/L;
    % MSE(i) = (S_normalized*w_t - S_normalized_t*w_t)'*(S_normalized*w_t - S_normalized_t*w_t)/L;
    % difference_scaled = (abs(S*w_t)/max(abs(S*w_t)) - abs(S*w)/max(abs(S*w)));
    % MSE_difference(i) = sum(difference_scaled)/L;
end

[~,index_min] = min(cost);
index_th = sorted_indices(index_min);
threshold = abs(w(index_th));
w_result = W(:,index_min);

% index_0 = sorted_indices(1:index_min);
% w_result = w;
% w_result(index_0) = 0;

end

% %% Iterative hard thresholding pursuit function
% function [cost,threshold,w_result,W,MSE] = iterativethreshold2(w,S,lambda)
% 
% N = length(w);
% cost = zeros(N,1);
% MSE = zeros(N,1);
% w_t = w;
% [~, sorted_indices] = sort(abs(w),'ascend');
% L = size(S,1);
% 
% W = zeros(N);
% 
% for i = 1:1:N
%     index = sorted_indices(i);
%     w_t(index) = 0;
%     % cost function evaluation
%     cost(i) = lambda*nnz(w_t)/N + (w - w_t)'*S'*S*(w - w_t)/L;
%     W(:,i) = w_t;
%     MSE(i) = (w - w_t)'*S'*S*(w - w_t)/L;
%     MSE_
% end
% 
% [~,index_min] = min(cost);
% index_th = sorted_indices(index_min);
% threshold = abs(w(index_th));
% w_result = W(:,index_min);
% 
% % index_0 = sorted_indices(1:index_min);
% % w_result = w;
% % w_result(index_0) = 0;
% 
% end