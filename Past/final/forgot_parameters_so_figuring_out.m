%% test what is the parameters used for the result shown in midterm

% initialize all parameters

%% set two desired patterns, one with beamwidth of 0.8, one with beamwidth of 0.08
left_win = -0.4;
right_win = 0.4;
desired_pattern08 = upperboundgen(win.leftend,win.rightend,[left_win right_win],L);       % The desired pattern
figure
plot(angle,desired_pattern08)
xlabel('Angle (degrees)');
ylabel('Amplitude')
title('Ideal beam pattern')

left_win = -0.04;
right_win = 0.04;
desired_pattern008 = upperboundgen(win.leftend,win.rightend,[left_win right_win],L);       % The desired pattern
figure
plot(angle,desired_pattern008)
xlabel('Angle (degrees)');
ylabel('Amplitude')
title('Ideal beam pattern')

error08 = norm(S*w - desired_pattern08,2)       % 4.6898e+03
error008 = norm(S*w - desired_pattern008,2)     % 1.4129e+03
% The generated pattern is more close to the reference pattern with
% beamwidth of 0.08. Draw the conclusion that the beamwidth was set to 0.08
% epsilon was set to 

%% modify the S matrix and reference pattern with loose range of 0.02 and 0.2
loose_range = 0.2;
desired_pattern_modified_loose02 = adjustRefpattern(desired_pattern008,[win.leftend win.rightend],[left_win right_win],L,loose_range);
S_new_loose02 = adjustSmatrix(S,[win.leftend win.rightend],[left_win right_win],L,loose_range);

loose_range = 0.02;
desired_pattern_modified_loose002 = adjustRefpattern(desired_pattern008,[win.leftend win.rightend],[left_win right_win],L,loose_range);
S_new_loose002 = adjustSmatrix(S,[win.leftend win.rightend],[left_win right_win],L,loose_range);

error008_loose02 = norm(S_new_loose02*w - desired_pattern_modified_loose02,2)       % 4.6898e+03
error008_loose002 = norm(S_new_loose002*w - desired_pattern_modified_loose002,2)     % 1.4129e+03
%%
loose_range = 0.03
desired_pattern_modified_loose = adjustRefpattern(desired_pattern008,[win.leftend win.rightend],[left_win right_win],L,loose_range);
S_new_loose = adjustSmatrix(S,[win.leftend win.rightend],[left_win right_win],L,loose_range);
error008_loose = norm(S_new_loose*w - desired_pattern_modified_loose,2)       % 4.6898e+03
