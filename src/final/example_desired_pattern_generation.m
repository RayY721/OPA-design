% Define the range for theta
theta = -18:0.01:18;

% Define the function
% The function is 1 when the absolute value of theta is less than 1, and 0 otherwise
f_theta = abs(theta) < 1;

% Plot the function
plot(theta, f_theta);
% grid on;
xlabel('\theta','fontsize',14);
ylabel('E(\theta)','fontsize',14);
title('The desired beampattern','fontsize',14);
axis([-18 18 -0.3 1.2]);