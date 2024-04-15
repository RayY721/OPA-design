function [W,E] = threshold(w)

N = length(w);
w_t = w;
[~, sorted_indices] = sort(abs(w),'ascend');
% L = size(S,1);
W = zeros(N);
E = zeros(N);
% S_app = diag(ones(N,1));
W(:,1) = w_t;
E(:,1) = w - w_t;
for i = 1:1:N-1
    index = sorted_indices(i);
    w_t(index) = 0;
    W(:,i+1) = w_t;
    E(:,i+1) = w - w_t;
end
% index_th = sorted_indices(index_min);
% threshold = abs(w(index_th));
end

% A question arised, how to tell the range of sidelobes, try something
% please
% maybe use the [-1,1] as a mainlobe region?

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