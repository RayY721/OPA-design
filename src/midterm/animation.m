%% Animation
% Create a figure
fig = figure();

% Define the range for the animation
numFrames = 1000;
theta = linspace(0, 2*pi, numFrames);

% Loop through each frame and update the plot
for i = 1:5:numFrames
    % Clear the previous plot
    clf
    
    % Create your animation content for this frame
    % x = cos(theta(i));
    % y = sin(theta(i));
    
    % Update the plot
    subplot(2,1,1)
    plot(angle,20*log10(abs(S*W(:,i))))
    yline(10*log10(MSE_app(i)))
    xlim([-18, 18]);
    ylim([-60, 0]);
    title("Array pattern with k =" + nnze(i));
    % xlabel('X-axis');
    % ylabel('Y-axis');

    % subplot(4,1,2)
    % plot(flip(MSE_new))
    % set(gca,'xdir','reverse')
    % xline(nnze(i))
    % title('MSE exact value')
    % text(nnze(i),MSE_new(i),"\leftarrow MSE = " +MSE_new(i))
    % 
    % subplot(4,1,3)
    % plot(flip(MSE_app))
    % set(gca,'xdir','reverse')
    % xline(nnze(i))
    % title('MSE approximation')
    % text(nnze(i),MSE_app(i),"\leftarrow MSE = " +MSE_app(i))

    subplot(2,1,2)
    plot(k,flip(nnze));
    set(gca,'xdir','reverse')
    xline(nnze(i))
    
    % Pause for a short time to control animation speed
    pause(0.01);
    
    % Capture the frame for the animation
    frame = getframe(fig);
    
    % Convert the frame to an image
    image = frame2im(frame);
    
    % Write the image to a video file (requires VideoWriter)
    if i == 1
        videoFileName = 'animation_with_enhanced.avi';
        vidObj = VideoWriter(videoFileName);
        open(vidObj);
    end
    writeVideo(vidObj, image);
end

% Close the video file
close(vidObj);

disp('Animation created successfully.');