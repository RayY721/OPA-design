% Define parameters
theta = -90:1:90;  % Angle range from -90 to 90 degrees
theta = linspace(-90,90,505);
mainlobe_width = 20; % Width of the main lobe in degrees
grating_lobe_width = 20; % Width of the grating lobes in degrees
% 
% % Compute the beam pattern
% mainlobe = sinc((theta / mainlobe_width) * pi); % Main lobe
% grating_lobe1 = sinc(((theta - 45) / grating_lobe_width) * pi); % Grating lobe 1
% grating_lobe2 = sinc(((theta + 45) / grating_lobe_width) * pi); % Grating lobe 2
% 
% % Combine the main lobe and grating lobes
% beam_pattern = mainlobe + grating_lobe1 + grating_lobe2;


% Normalize the beam pattern
beam_pattern = beam_pattern / max(beam_pattern);

% Plot the beam pattern
figure;
plot(theta, beam_pattern, 'LineWidth', 2);
xlabel('Angle (degrees)');
ylabel('Amplitude');
title('Beam pattern');
grid on;
axis([-90 90 0 1]);