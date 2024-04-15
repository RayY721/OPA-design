%% Plot the Gaussian Curve
gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
xx = linspace(-5,5,1000);
mu = 0;
sig = 2.3;
amp = 0.006;
vo = 0; 
y = gaus(xx,mu,sig,amp,vo);
% Plot gaussian
figure
stem(abs(W_rew2(:,5)))
hold on
plot(y, 'b-', 'LineWidth',3)