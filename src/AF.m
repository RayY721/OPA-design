% This function takes inputs as:
% elementPos - position of elements
% leftend
% rightend
% Res - the # of samples within the range defined by leftend and rightend
% k - wavenumber, spatial frequency
% ang - target angle
% plot_index - wether to plot. (true to plot, false not to plot)
% w - weight of each element(only considering the logic 0/1 weight now)
function AF_dB = AF(elementPos,leftend,rightend,Res,k,ang,plot_index,w)
if nargin < 8
    w = ones(size(elementPos));
end
activeelementPos = elementPos(w~=0);
N = numel(activeelementPos);
% ang - The desired angle, a scaler in 1D case,
% deviate from the vertical, positive steer right, negative steer left
angle = linspace(leftend,rightend,Res);
% angle = linspace(-90,90,Res);
theta = angle*pi/180;
ang_rad = ang*pi/180;
% u = [sin(ang),cos(ang)];
% The desired phase shift on mth element
D_phi = k*sin(ang_rad)*activeelementPos;   
% The simplest beamformer for the desired direction
w = exp(1i*D_phi);      % given the desired angle, the delay between origin and the mth element

% The array response
% exp(1i*k*cos(theta(i))*elementPos)
S = exp(1i*k*activeelementPos*sin(theta))';
AF = S*w/N;
% AF = S*w;
AF_dB = 20*log10(abs(AF));          % AF suppose to be magnitude
% AF_dB = pow2db(abs(AF));
sine_angle = sin(theta);
if plot_index == true
    figure
    plot(sine_angle,AF_dB)
    % ylim([-40 0])
    hold on 
    % line([ang ang], ylim, 'Color', 'r', 'LineStyle', '--');
end
end