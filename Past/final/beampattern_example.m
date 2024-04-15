%% Beampattern example
close all
clc
clear
% Parameters
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;
N = 11;                % to have a aperture of 1020 lambda
d = 1.5*param.lambda;
L = 5000;             % The resolution is at least 0.035 degree
fov.left = -90;
fov.right = 90;
Res = (fov.right - fov.left)/(L - 1);
angle = linspace(fov.left,fov.right,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix
w = ones(N,1);

figure;
% plot(angle,20*log10(abs(S*w)/max(abs(S*w))),'LineWidth', 2)
polarplot(theta,abs(S*w)/max(abs(S*w)),'LineWidth', 2)
% Assuming you've already created your polar plot using the polarplot function
pax = gca;  % Get current polar axes
pax.ThetaZeroLocation = 'top';
% pax.ThetaDirection = 'clockwise';
pax.ThetaLim = [-90 90];  % Limit the display to the upper hemisphere
pax.ThetaAxis.FontSize = 12;  % Set the font size for the angular axis
pax.RAxis.FontSize = 12;



%% Beampattern example (normal plot)
close all
clc
clear
% Parameters
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;
N = 11;                % to have a aperture of 1020 lambda

d = 1.5*param.lambda;
L = 5000;             % The resolution is at least 0.035 degree
fov.left = -90;
fov.right = 90;
Res = (fov.right - fov.left)/(L - 1);
angle = linspace(fov.left,fov.right,L)';  % in degree
theta = angle*pi/180;
x = (-(N-1)*d/2:d:(N-1)*d/2)';
S = exp(1i*param.k*sin(theta)*x');            % S matrix

w = ones(N,1);
% beam_pattern = beam_pattern / max(beam_pattern);

figure;
% plot(angle,20*log10(abs(S*w)))
plot(angle,20*log10(abs(S*w)/max(abs(S*w))),'LineWidth', 2);
xlabel('Angle (degrees)','fontsize',14);
ylabel('Intensity(dB)','fontsize',14);
title('Beampattern','fontsize',14);
% grid on;
xlim([-90 90])
ylim([-70 0])


%% Beampattern example (gradient of N and spacing)
close all
clc
clear
% Parameters
param.c = physconst('LightSpeed');
param.lambda = 1550e-9;
param.fc = param.c/param.lambda;
param.k = 2*pi/param.lambda;
N1 = 11;                % to have a aperture of 1020 lambda
N2 = 36;
N3 = 51;
d1 = 0.5*param.lambda;
d2 = 1.5*param.lambda;
d3 = 2.5*param.lambda;
L = 5000;             % The resolution is at least 0.035 degree
fov.left = -90;
fov.right = 90;
Res = (fov.right - fov.left)/(L - 1);
angle = linspace(fov.left,fov.right,L)';  % in degree
theta = angle*pi/180;
x11 = (-(N1-1)*d1/2:d1:(N1-1)*d1/2)';
x21 = (-(N2-1)*d1/2:d1:(N2-1)*d1/2)';
x31 = (-(N3-1)*d1/2:d1:(N3-1)*d1/2)';
% x = (-(N-1)*d/2:d:(N-1)*d/2)';
x12 = (-(N1-1)*d2/2:d2:(N1-1)*d2/2)';
x13 = (-(N1-1)*d3/2:d3:(N1-1)*d3/2)';

S11 = exp(1i*param.k*sin(theta)*x11');            % S matrix
S21 = exp(1i*param.k*sin(theta)*x21');            % S matrix
S31 = exp(1i*param.k*sin(theta)*x31');            % S matrix
% S11 = exp(1i*param.k*sin(theta)*x11');            % S matrix
S12 = exp(1i*param.k*sin(theta)*x12');            % S matrix
S13 = exp(1i*param.k*sin(theta)*x13');            % S matrix

w1 = ones(N1,1);
w2 = ones(N2,1);
w3 = ones(N3,1);
% beam_pattern = beam_pattern / max(beam_pattern);

figure;
% plot(angle,20*log10(abs(S*w)))
a1 = plot(angle,20*log10(abs(S11*w1)/max(abs(S11*w1))),'LineWidth', 2);
M1 = "N = 11";
hold on 
a2 = plot(angle,20*log10(abs(S21*w2)/max(abs(S21*w2))),'LineWidth', 2);
M2 = "N = 31";
a3 = plot(angle,20*log10(abs(S31*w3)/max(abs(S31*w3))),'LineWidth', 2);
M3 = "N = 51";
set(gca, 'FontSize', 24);
xlabel('\theta','fontsize',37);
ylabel('Intensity(dB)','fontsize',36);
title('Beampattern with half-wavelength spacing','fontsize',36);
% grid on;
legend([a1,a2,a3],[M1,M2,M3],'FontSize',39)
xlim([-90 90])
ylim([-50 0])

figure;
% plot(angle,20*log10(abs(S*w)))
a1 = plot(angle,20*log10(abs(S11*w1)/max(abs(S11*w1))),'LineWidth', 2);
M1 = "d = 0.5 wavelength";
hold on 
a2 = plot(angle,20*log10(abs(S12*w1)/max(abs(S12*w1))),'LineWidth', 2);
M2 = "d = 1.5 wavelength";
a3 = plot(angle,20*log10(abs(S13*w1)/max(abs(S13*w1))),'LineWidth', 2);
M3 = "d = 2.5 wavelength";
set(gca, 'FontSize', 24);
xlabel('\theta','fontsize',37);
ylabel('Intensity(dB)','fontsize',36);
title('Beampattern with N = 11','fontsize',36);
% grid on;
legend([a1,a2,a3],[M1,M2,M3],'FontSize',39)
xlim([-90 90])
ylim([-50 0])
