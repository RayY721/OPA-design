%% this script make the plots for the thesis. The excitation is generated with epsilon of 0.06,
% transition region of 0.02
figure
subplot(3,1,1)
plot(angle,20*log10(abs(S*w)))
subplot(3,1,2)
stem(abs(w))
title('The amplitude distribution')
ylabel('Amplitude')
xlabel('Index')
subplot(3,1,3)
semilogy(flip(sort(abs(w))))
title('The sorted amplitude distribution')
ylabel('Amplitude(dB)')
xlabel('Index')

%%
figure;
plot(angle,20*log10(abs(S_normalized*W_rew(:,1))))
hold on
plot(angle,20*log10(abs(S_normalized*W_rew(:,2))))
plot(angle,20*log10(abs(S_normalized*W_rew(:,5))))