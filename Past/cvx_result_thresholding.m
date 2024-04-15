load('complex_weight.mat', 'W_modified_small')

w = W_modified_small(:,4);
figure
subplot(2,2,1)
plot(abs(w));
title('Original weights')

subplot(2,2,2)
indexn = abs(w) <= 0.01;
% w(indexn == 0,4);
w_thresheld = w;
w_thresheld(indexn == 1) = 0;
plot(abs(w_thresheld))
title('Thresheld weights')

subplot(2,2,3)
plot(angle,20*log10(abs(S*w)))
title('AF before thresholding')

subplot(2,2,4)
plot(angle,20*log10(abs(S*w_thresheld)))
title('AF after thresholding')

norm(abs(S*w_thresheld - S*w))
% norm(abs(S*(w_thresheld-w)))
% norm(S*(w_thresheld-w))
figure
semilogy(sort(abs(w)))

%% check the structure of the S matrix
figure
mesh(abs(S'*S))
figure
mesh(abs(S*S'))