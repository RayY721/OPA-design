%% Iterative hard thresholding pursuit function
function [cost,threshold,w_result,W,MSE] = iterativethreshold(w,S,lambda)

N = length(w);
cost = zeros(N,1);
MSE = zeros(N,1);
w_t = w;
[~, sorted_indices] = sort(abs(w),'ascend');
L = size(S,1);

W = zeros(N);

for i = 1:1:N
    index = sorted_indices(i);
    w_t(index) = 0;
    % cost function evaluation
    cost(i) = lambda*nnz(w_t)/N + (w - w_t)'*S'*S*(w - w_t)/L;
    W(:,i) = w_t;
    MSE(i) = (w - w_t)'*S'*S*(w - w_t)/L;
end

[~,index_min] = min(cost);
index_th = sorted_indices(index_min);
threshold = abs(w(index_th));
w_result = W(:,index_min);

% index_0 = sorted_indices(1:index_min);
% w_result = w;
% w_result(index_0) = 0;

end