% Load two frames (images) from your video or dataset
% For example, frame50 and frame51
v = VideoReader("DATA/closedcellday_2022_09_06.mp4");
video = read(v);
red = squeeze(video(:,:,1,:)); % Extract the red channel

% Full-res Frames 3 and 4. Flip arrays, images are upside down
fframe = red(1:1000,600:1450,3) ; % for real data image, 3rd image in video
fframe2= red(1:1000,600:1450,4) ; DT = 30*60; % 4th image in video, 1/2 h later
% coordinate arrays at full resolution
fy = 1:size(fframe,1);
fx = 1:size(fframe,2);

% resize images (arrays) to make computations faster  
% shrinkfactor=1;
shrinkfactor=5; invshrinkfactor = 0.2;
frame = imresize(fframe , invshrinkfactor);
frame2= imresize(fframe2, invshrinkfactor);


% Ensure the frames are in double precision
frame = double(frame);
frame2 = double(frame2);

%%% Wavelet transform inputs (scales and angles) 
% 24 Angles from 0 to pi 
NSCALES = 24;
Angles = 0:pi/NSCALES:pi ;

% a LOGARITHMIC set of 10 Scales, shrinkfactor resizes them
Scales = 10.^(1:.05:1.9);
Scales = Scales/shrinkfactor; 


% Perform 2D Continuous Wavelet Transform on both frames
cwtCauchy1 = cwtft2(frame, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
spec1 = squeeze(cwtCauchy1.cfs);

cwtCauchy2 = cwtft2(frame2, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
spec2 = squeeze(cwtCauchy2.cfs);


% Compute the cross-wavelet spectrum (XWT)
xwt = spec1 .* conj(spec2);

% Compute coherence (amplitude) and phase difference
coherence = abs(xwt);
phase_difference = angle(xwt);

clear('video','xwt','cwtCauchy1','cwtCauchy2')
%%
plot_waverose_Power(spec1,Scales)
plot_waverose_Power(spec2,Scales)


%%
peak_list = plot_waverose_Coherence_Maximums(coherence,Scales);

%%
%{
% Initialize a figure to plot multiple coherence images
figure;

% Define common fractions of pi for displaying angles
pi_fractions = {'0', '\pi/6', '\pi/4', '\pi/3', '\pi/2', '2\pi/3', '3\pi/4', '5\pi/6', '\pi', ...
                '7\pi/6', '5\pi/4', '4\pi/3', '3\pi/2', '5\pi/3', '7\pi/4', '11\pi/6', '2\pi'};
pi_fraction_values = linspace(0, 2*pi, length(pi_fractions));  % Corresponding numerical values for ticks

% Loop through each pair of scale and angle in the peak_list
for i = 1:size(peak_list, 2)
    
    % Extract the real scale and angle from the peak_list
    real_scale = peak_list(1, i);  % Scale from peak_list
    real_angle = peak_list(2, i);  % Angle from peak_list
    
    % Find the index of the closest scale in the Scales array
    [~, scale_idx] = min(abs(Scales - real_scale));
    
    % Define an Angles array corresponding to the indices in the coherence 4th dimension
    % Assume Angles_full is a vector of angle values (e.g., linspace(0, 2*pi, num_angles))
    num_angles = size(coherence, 4);  % Number of angles in coherence
    Angles = linspace(0, 2*pi, num_angles);  % Create the angle array based on your data structure
    
    % Find the index of the closest angle in the Angles array
    [~, angle_idx] = min(abs(Angles - real_angle));
    
    % Extract the corresponding coherence slice for the current scale and angle
    coherence_slice = coherence(:, :, scale_idx, angle_idx);
    
    % Convert the real_angle to a fraction of pi for the title
    [~, angle_fraction_idx] = min(abs(pi_fraction_values - real_angle));  % Find closest pi fraction
    angle_str = pi_fractions{angle_fraction_idx};  % Get the corresponding fraction of pi string
    
    % Plot the coherence using imagesc
    subplot(ceil(sqrt(size(peak_list, 2))), ceil(sqrt(size(peak_list, 2))), i);  % Subplot for multiple plots
    imagesc(coherence_slice);
    
    % Set title with scale and angle in fractions of pi
    title(['Scale: ', num2str(real_scale), ', Angle: ', angle_str]);
    
    % Customize the colorbar to display ticks from -pi to pi
    c = colorbar;
    
    % Set axis equal for consistent plotting
    axis equal;
end

% Customize the overall figure title
sgtitle('Coherence for Each Scale/Angle Peak');

%%
% Initialize a figure to plot multiple phase_difference images
figure;

% Define common fractions of pi for displaying angles
pi_fractions = {'0', '\pi/6', '\pi/4', '\pi/3', '\pi/2', '2\pi/3', '3\pi/4', '5\pi/6', '\pi', ...
                '7\pi/6', '5\pi/4', '4\pi/3', '3\pi/2', '5\pi/3', '7\pi/4', '11\pi/6', '2\pi'};
pi_fraction_values = linspace(0, 2*pi, length(pi_fractions));  % Corresponding numerical values for ticks

% Loop through each pair of scale and angle in the peak_list
for i = 1:size(peak_list, 2)
    
    % Extract the real scale and angle from the peak_list
    real_scale = peak_list(1, i);  % Scale from peak_list
    real_angle = peak_list(2, i);  % Angle from peak_list
    
    % Find the index of the closest scale in the Scales array
    [~, scale_idx] = min(abs(Scales - real_scale));
    
    % Define an Angles array corresponding to the indices in the phase_difference 4th dimension
    % Assume Angles_full is a vector of angle values (e.g., linspace(0, 2*pi, num_angles))
    num_angles = size(phase_difference, 4);  % Number of angles in phase_difference
    Angles = linspace(0, 2*pi, num_angles);  % Create the angle array based on your data structure
    
    % Find the index of the closest angle in the Angles array
    [~, angle_idx] = min(abs(Angles - real_angle));
    
    % Extract the corresponding phase_difference slice for the current scale and angle
    phase_slice = phase_difference(:, :, scale_idx, angle_idx);
    
    % Convert the real_angle to a fraction of pi for the title
    [~, angle_fraction_idx] = min(abs(pi_fraction_values - real_angle));  % Find closest pi fraction
    angle_str = pi_fractions{angle_fraction_idx};  % Get the corresponding fraction of pi string
    
    % Plot the phase_difference using imagesc
    subplot(ceil(sqrt(size(peak_list, 2))), ceil(sqrt(size(peak_list, 2))), i);  % Subplot for multiple plots
    imagesc(phase_slice);
    
    % Set title with scale and angle in fractions of pi
    title(['Scale: ', num2str(real_scale), ', Angle: ', angle_str]);
    
    % Customize the colorbar to display ticks from -pi to pi
    c = colorbar;
    clim([-pi pi]);  % Set color axis limits from -pi to pi
    set(c, 'YTick', -pi:pi/2:pi, 'YTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});  % Set ticks as fractions of pi
    
    % Set axis equal for consistent plotting
    axis equal;
end

% Customize the overall figure title
sgtitle('Phase Difference for Each Scale/Angle Peak');

%}

%% Part 1: Plot coherence with red contour based on the 60% threshold
figure;

% Define common fractions of pi for displaying angles
pi_fractions = {'0', '\pi/6', '\pi/4', '\pi/3', '\pi/2', '2\pi/3', '3\pi/4', '5\pi/6', '\pi', ...
                '7\pi/6', '5\pi/4', '4\pi/3', '3\pi/2', '5\pi/3', '7\pi/4', '11\pi/6', '2\pi'};
pi_fraction_values = linspace(0, 2*pi, length(pi_fractions));  % Corresponding numerical values for ticks

% Preallocate a third row for the phase mean values in peak_list
peak_list(3, :) = NaN;

% Preallocate a cell array to store coherence masks
coherence_masks = cell(1, size(peak_list, 2));

for i = 1:size(peak_list, 2)
    
    % Extract the real scale and angle from the peak_list
    real_scale = peak_list(1, i);  % Scale from peak_list
    real_angle = peak_list(2, i);  % Angle from peak_list
    
    % Find the index of the closest scale in the Scales array
    [~, scale_idx] = min(abs(Scales - real_scale));
    
    % Define an Angles array corresponding to the indices in the coherence 4th dimension
    num_angles = size(coherence, 4);  % Number of angles in coherence
    Angles = linspace(0, 2*pi, num_angles);  % Create the angle array based on your data structure
    
    % Find the index of the closest angle in the Angles array
    [~, angle_idx] = min(abs(Angles - real_angle));
    
    % Extract the corresponding coherence slice for the current scale and angle
    coherence_slice = coherence(:, :, scale_idx, angle_idx);
    
    % Define the threshold for the top 90% of the max coherence value
    coherence_max = max(coherence_slice(:));
    coherence_threshold = 0.6 * coherence_max;
    
    % Create the mask where coherence is above the threshold
    coherence_mask = coherence_slice >= coherence_threshold;
    
    % Store the mask in the cell array
    coherence_masks{i} = coherence_mask;  % Store mask for later use
    
    % Plot the coherence using imagesc
    subplot(ceil(sqrt(size(peak_list, 2))), ceil(sqrt(size(peak_list, 2))), i);  % Subplot for multiple plots
    imagesc(coherence_slice);
    hold on;
    
    % Overlay the red contour for the coherence mask
    [contour_x, contour_y] = find(coherence_mask);
    scatter(contour_y, contour_x, 'r', 'LineWidth', 1);  % Contour plot in red
    
    % Convert the real_angle to a fraction of pi for the title
    [~, angle_fraction_idx] = min(abs(pi_fraction_values - real_angle));  % Find closest pi fraction
    angle_str = pi_fractions{angle_fraction_idx};  % Get the corresponding fraction of pi string
    
    % Set title with scale and angle in fractions of pi
    title(['Scale: ', num2str(real_scale), ', Angle: ', angle_str]);
    
    % Customize the colorbar
    colorbar;
    
    % Set axis equal for consistent plotting
    axis equal;
end

sgtitle('Coherence with Mask for Each Scale/Angle Peak');

%% Part 2: Plot phase difference with masked values and compute mean
figure;

for i = 1:size(peak_list, 2)
    
    % Extract the real scale and angle from the peak_list
    real_scale = peak_list(1, i);  % Scale from peak_list
    real_angle = peak_list(2, i);  % Angle from peak_list
    
    % Find the index of the closest scale in the Scales array
    [~, scale_idx] = min(abs(Scales - real_scale));
    
    % Define an Angles array corresponding to the indices in the phase_difference 4th dimension
    num_angles = size(phase_difference, 4);  % Number of angles in phase_difference
    Angles = linspace(0, 2*pi, num_angles);  % Create the angle array based on your data structure
    
    % Find the index of the closest angle in the Angles array
    [~, angle_idx] = min(abs(Angles - real_angle));
    
    % Extract the corresponding phase_difference slice for the current scale and angle
    phase_slice = phase_difference(:, :, scale_idx, angle_idx);
    
    % Retrieve the corresponding mask from Part 1
    coherence_mask = coherence_masks{i};  % Get the correct mask for this scale/angle
    
    % Apply the mask from coherence to the phase slice (mask where phase should be visible)
    phase_masked = phase_slice;
    phase_masked(~coherence_mask) = NaN;  % Set non-coherent areas to NaN
    
    % Plot the phase_difference using imagesc
    subplot(ceil(sqrt(size(peak_list, 2))), ceil(sqrt(size(peak_list, 2))), i);  % Subplot for multiple plots
    imagesc(phase_masked);
    
    % Set title with scale and angle in fractions of pi
    title(['Scale: ', num2str(real_scale), ', Angle: ', angle_str]);
    
    % Customize the colorbar to display ticks from -pi to pi
    c = colorbar;
    clim([-pi pi]);  % Set color axis limits from -pi to pi
    set(c, 'YTick', -pi:pi/2:pi, 'YTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});  % Set ticks as fractions of pi
    
    % Set axis equal for consistent plotting
    axis equal;
    
    % Calculate the mean phase value for the masked region
    phase_mean = mean(phase_slice(coherence_mask), 'omitnan');
    
    % Add the phase mean value to the third row of peak_list
    peak_list(3, i) = phase_mean;
end

% Customize the overall figure title
sgtitle('Phase Difference with Masked Regions for Each Scale/Angle Peak');

% Display the updated peak_list
disp('Updated peak_list with mean phase values:');
disp(peak_list);


%% Calculate Speed from Phase Shift and Plot Sine Waves

% Constants
pixel_size_km = 2 * shrinkfactor;       % Each pixel equals 2 km
time_interval_s = 1800;  % 30 minutes in seconds

% Initialize arrays to store results
num_peaks = size(peak_list, 2);
speeds = zeros(1, num_peaks);
distance_shifts = zeros(1, num_peaks);

for i = 1:num_peaks
    % Extract scale and mean phase difference from peak_list
    scale = peak_list(1, i);
    mean_phase_difference = peak_list(3, i);
    
    % Calculate the wavelength in km
    % Since scale is half the wavelength in pixels, and each pixel is 2 km
    wavelength_km = scale * 2 * pixel_size_km;
    
    % Calculate the distance shift in km
    distance_shift_km = mean_phase_difference * wavelength_km / (2 * pi);
    
    % Calculate the speed in km/s
    speed_km_per_s = distance_shift_km / time_interval_s;
    
    % Store the results
    distance_shifts(i) = distance_shift_km;
    speeds(i) = speed_km_per_s;
    
    % Plot the sine waves
    x_values = linspace(0, wavelength_km, 100);
    sine_wave_original = sin(2 * pi * x_values / wavelength_km);
    sine_wave_with_phase = sin(2 * pi * x_values / wavelength_km + mean_phase_difference);
    
    figure;
    plot(x_values, sine_wave_original, 'b-', 'LineWidth', 2);
    hold on;
    plot(x_values, sine_wave_with_phase, 'r--', 'LineWidth', 2);
    xlabel('Distance (km)');
    ylabel('Amplitude');
    legend('Original Sine Wave', 'Shifted Sine Wave');
    title(sprintf('Sine Wave with Phase Shift for Peak %d\nScale: %.2f pixels, Wavelength: %.2f km\nPhase: %.4f radians (%.4f\\pi)\nShift: %.4f km, Speed: %.4f m/s', ...
    i, scale, wavelength_km, mean_phase_difference, mean_phase_difference / pi, distance_shift_km, speed_km_per_s * 1000));
    grid on;
    xlim([0, wavelength_km]);
    hold off;
end

% Add the speeds to the peak_list
peak_list(4, :) = speeds*1000;

% Display the updated peak_list
disp('Updated peak_list with mean phase values and speeds:');
disp('Rows: 1-Scale, 2-Angle, 3-Mean Phase Difference, 4-Speed (km/s)');
disp(peak_list);


%% Plot Rose Plot with Power and Coherence + pikes

% Assuming 'spec1' and 'coherence' data is already computed in your workspace.
% Compute the wavelet power
power = abs(spec1).^2;

% Define buffer to avoid edge effects
buffer = round(max(Scales));

% Calculate the mean inner power within the buffered region
innerpower = squeeze(mean(mean(power(buffer:size(power, 1) - buffer, ...
                                      buffer:size(power, 2) - buffer, :, :))));

% Calculate the mean of innerpower across scales
meanbyscale = squeeze(mean(transpose(innerpower)));

% Initialize the angle spectrum container (same size as innerpower)
anglespec_power = innerpower .* 0;
for isc = 1:24   % Loop over angles (assuming 24 angular bins)
    anglespec_power(:, isc) = squeeze(innerpower(:, isc)) ./ transpose(meanbyscale);
end

%% Calculate and Plot Coherence (for lower half)

% Calculate inner coherence using the same buffer
innercoherence = squeeze(mean(mean(coherence(buffer:size(coherence, 1) - buffer, ...
                                            buffer:size(coherence, 2) - buffer, :, :))));

% Calculate the mean of innercoherence across scales
mean_coherence_byscale = squeeze(mean(transpose(innercoherence)));

% Initialize the angle spectrum container for coherence
anglespec_coherence = innercoherence .* 0;
for isc = 1:24   % Loop over angles
    anglespec_coherence(:, isc) = squeeze(innercoherence(:, isc)) ./ transpose(mean_coherence_byscale);
end

%% Plot the Upper Half (Power Spectrum) (0 to pi)
% Define angles for the upper half (from 0 to pi)
Angles_pos = linspace(0, pi, size(anglespec_power, 2));

% Create meshgrid for polar coordinates (upper half)
[Theta_pos, R_pos] = meshgrid(Angles_pos, Scales);

% Convert polar coordinates to Cartesian for plotting with pcolor (upper half)
[X_pos, Y_pos] = pol2cart(Theta_pos, R_pos);

% Create a figure and first axis for the upper half
figure;
ax1 = axes;
pcolor(ax1, X_pos, Y_pos, anglespec_power);  % Plot power spectrum (upper half)
shading interp;  % Smooth shading for the plot
colormap(ax1);  % Use the 'parula' colormap for the upper half
%colorbar(ax1, 'eastoutside');  % Add colorbar for the first plot
title('Upper Half Wavelet Power Rose Plot');
axis equal;  % Ensure equal scaling
set(ax1, 'Position', [0.1, 0.1, 0.75, 0.75]);  % Set axis position
ax1.XTick = [];  % Remove x-axis ticks
ax1.YTick = [];  % Remove y-axis ticks

% Define angles for the lower half (from pi to 2*pi)
Angles_neg = linspace(pi, 2*pi, size(anglespec_coherence, 2));

% Create meshgrid for polar coordinates (lower half)
[Theta_neg, R_neg] = meshgrid(Angles_neg, Scales);

% Convert polar coordinates to Cartesian for plotting with pcolor (lower half)
[X_neg, Y_neg] = pol2cart(Theta_neg, R_neg);

% Create second axis for the lower half and plot it
ax2 = axes;
pcolor(ax2, X_neg, Y_neg, anglespec_coherence);  % Plot coherence spectrum (lower half)
shading interp;  % Smooth shading
colormap(ax2, 'autumn');  % Use a different colormap ('jet') for the lower half
linkaxes([ax1,ax2])
axis equal;  % Ensure equal scaling
set(ax2, 'Position', [0.1, 0.1, 0.75, 0.75]);  % Ensure ax2 has the same position as ax1
set(ax2, 'Color', 'none');  % Make the second axis transparent
ax2.XTick = [];  % Remove x-axis ticks
ax2.YTick = [];  % Remove y-axis ticks

xlim([min(X_pos,[],"all")-1 max(X_pos,[],"all")+1 ])
ylim([-max(Y_pos,[],"all")-1 max(Y_pos,[],"all")+1 ])


%% Overlay Peaks on Lower Half (if peak_list is defined)

max_scales = peak_list(1, :);  % Scale values from peak_list
max_angles = peak_list(2, :);  % Angle values from peak_list
speeds = peak_list(4, :);      % Speed values from peak_list

% Convert peaks' polar coordinates to Cartesian coordinates for lower half
[max_X, max_Y] = pol2cart(max_angles+pi, max_scales);

% Overlay the peaks on the lower half plot
hold(ax2, 'on');
hold(ax1, 'on');
plot(ax2, max_X, max_Y, 'k*', 'MarkerSize', 10);  % Plot peaks as black stars

% Annotate each peak with its speed value (in m/s)
for i = 1:length(max_X)
    text(ax2, max_X(i) * 1.05, max_Y(i) * 1.05, sprintf('%.2f m/s', speeds(i)), 'Color', 'k', 'FontSize', 10);
end


%% Add Radial Rings Using Actual Scales

% Add radial rings to both halves, using actual scales
ring_radii = Scales;  % Use scales as radii for the rings

for i = 1:length(ring_radii)
    if ~ismember(i, [2:3, 5:10, 12:13, 15, 17])  % Skip specific rings for clarity
        theta_ring = linspace(0, 2*pi, 100);  % Full circle for each ring
        [x_ring, y_ring] = pol2cart(theta_ring, ring_radii(i));  % Convert to Cartesian coordinates
        
        % Plot rings on both axes (ax1 for upper half, ax2 for lower half)
        plot(ax1, x_ring, y_ring, 'k--');  % Plot the radial ring on the upper half
        plot(ax2, x_ring, y_ring, 'k--');  % Plot the radial ring on the lower half
        
        % Add text for scale labels on the right side
        text(max(x_ring) * 1.05, 0, num2str(ring_radii(i), '%.2f'), 'HorizontalAlignment', 'left');
    end
end

%% Display Angle Ticks for Both Halves

% Upper half angular ticks
angle_ticks_pos = linspace(0, pi, 7);  % Ticks for the upper half
angle_labels_pos = {'0', '\pi/6', '\pi/3', '\pi/2', '2\pi/3', '5\pi/6', '\pi'};

for i = 1:length(angle_ticks_pos)
    angle_rad = angle_ticks_pos(i);
    plot(ax1, [0, max(X_pos(:)) * cos(angle_rad)], [0, max(Y_pos(:)) * sin(angle_rad)], 'k--');
    text(ax1, max(X_pos(:)) * cos(angle_rad) * 1.1, max(Y_pos(:)) * sin(angle_rad) * 1.1, ...
        angle_labels_pos{i}, 'HorizontalAlignment', 'center');
end

% Lower half angular ticks
angle_ticks_neg = linspace(pi, 2*pi, 7);  % Ticks for the lower half
angle_labels_neg = {'\pi', '7\pi/6', '4\pi/3', '3\pi/2', '5\pi/3', '11\pi/6', '2\pi'};
    
    % Display angle ticks for the lower half
    for i = 1:length(angle_ticks_neg)-1
        angle_rad = angle_ticks_neg(i);
        
        % Lower half
        plot([0 max(X_pos(:))*cos(angle_rad)], [0 max(Y_pos(:))*sin(angle_rad)], 'k--');
        text(max(X_pos(:))*cos(angle_rad)*1.1, max(Y_pos(:))*sin(angle_rad)*1.1, ...
            angle_labels_neg{i}, 'HorizontalAlignment', 'center');
    end
   
original_pos = get(ax1, 'Position');
c1 = colorbar(ax1, 'eastoutside');
c1_pos = get(c1, 'Position');  % Get the position of the colorbar
c1_pos(1) = c1_pos(1) + 0.05;  % Shift it further to the right (increase the x-position)
set(c1, 'Position', c1_pos);   % Set the new position
set(ax1, 'Position', original_pos);
ylabel(c1, 'Power');  

original_pos = get(ax1, 'Position');
c2 = colorbar(ax2, 'westoutside');
c2_pos = get(c2, 'Position');  % Get the position of the colorbar
c2_pos(1) = c2_pos(1) - 0.05;  % Shift it further to the right (increase the x-position)
set(c2, 'Position', c2_pos);   % Set the new position
set(ax2, 'Position', original_pos);
ylabel(c2, 'Cross-Wavelet Power');  


%%

function plot_waverose_Power(spec1,Scales)
    power = abs(spec1) .^2;
    buffer = round( max(Scales) );
    innerpower = squeeze( mean(mean( power(buffer:size(power,1)-buffer, ...
                                           buffer:size(power,2)-buffer, :,:) )));

    meanbyscale = squeeze( mean(transpose(innerpower)) );

    anglespec = innerpower .* 0;     % right sized container for angle spectrum
    for isc = 1:24   %size( transpose(Scales) )
        anglespec(:,isc) = squeeze(innerpower(:,isc)) ./ transpose(meanbyscale);
    end 
    
    % Define scales and angles for the upper half
    Angles_pos = linspace(0, pi, size(anglespec, 2));  % Angles from 0 to pi (upper half)
    % Scales = 10.^(1:.05:1.9);  % Logarithmic scales (adjust based on your data)
    
    % Extend angles to cover the full circle (0 to 2*pi)
    Angles_full = [Angles_pos, Angles_pos + pi];  % Full angles from 0 to 2*pi
    
    % Create meshgrid for polar coordinates
    [Theta, R] = meshgrid(Angles_full, Scales);  % Full angle grid for both halves
    
    % Extend the data to cover the lower half, not just mirrored
    anglespec_full = [anglespec, anglespec];  % Repeat the data for the second half (adjust based on physical symmetry)
    
    % Convert polar coordinates to Cartesian for plotting with pcolor
    [X, Y] = pol2cart(Theta, R);
    
    % Plot using Cartesian coordinates with pcolor
    figure;
    pcolor(X, Y, anglespec_full);  % Use pcolor to plot the wavelet power data
    shading interp;  % Smooth shading
    
    % Customize the plot with labels and title
    title('Wavelet Power Full Rose Plot');
    colorbar;  % Add colorbar to indicate wavelet power
    ylabel(colorbar, 'Wavelet Power');  % Add label to the colorbar
    
    % Remove Cartesian ticks 
    set(gca, 'XTick', [], 'YTick', []);  
    axis equal;  % Ensure equal scaling
    
    % Add angular lines with fractions of pi as labels
    hold on;
    angle_ticks_pos = linspace(0, pi, 7);  % Upper half ticks (0 to pi)
    angle_ticks_neg = linspace(pi, 2*pi, 7); % Lower half ticks (pi to 2*pi)
    angle_labels_pos = {'0', '\pi/6', '\pi/3', '\pi/2', '2\pi/3', '5\pi/6', '\pi'};  % Labels for upper half
    angle_labels_neg = {'\pi', '7\pi/6', '4\pi/3', '3\pi/2', '5\pi/3', '11\pi/6', '2\pi'};  % Labels for lower half
    
    % Display angle ticks for the upper half
    for i = 1:length(angle_ticks_pos)
        angle_rad = angle_ticks_pos(i);
        
        % Upper half
        plot([0 max(X(:))*cos(angle_rad)], [0 max(Y(:))*sin(angle_rad)], 'k--');
        text(max(X(:))*cos(angle_rad)*1.1, max(Y(:))*sin(angle_rad)*1.1, ...
            angle_labels_pos{i}, 'HorizontalAlignment', 'center');
    end
    
    % Display angle ticks for the lower half
    for i = 1:length(angle_ticks_neg)-1
        angle_rad = angle_ticks_neg(i);
        
        % Lower half
        plot([0 max(X(:))*cos(angle_rad)], [0 max(Y(:))*sin(angle_rad)], 'k--');
        text(max(X(:))*cos(angle_rad)*1.1, max(Y(:))*sin(angle_rad)*1.1, ...
            angle_labels_neg{i}, 'HorizontalAlignment', 'center');
    end
    
    % Add radial rings using actual Scales values
    ring_radii = Scales;  % Use actual scales for the rings
    for i = 1:length(ring_radii)
        if ~ismember(i, [2:3 ,5:10,12:13,15,17])
            theta_ring = linspace(0, 2*pi, 100);  % Full circle for each ring
            [x_ring, y_ring] = pol2cart(theta_ring, ring_radii(i));
            plot(x_ring, y_ring, 'k--');  % Plot the radial ring
    
            text(max(x_ring)*1.05, 0, num2str(ring_radii(i), '%.2f'), 'HorizontalAlignment', 'left')
        end
    end
    hold off;
end

function plot_waverose_Coherence(coherence,Scales)

    innercoherence = squeeze(  mean(mean( coherence(50:800,200:800, :,:) ))  );
    meanbyscale = squeeze( mean(transpose(innercoherence)) );

    
    anglespec = innercoherence .* 0;     % right sized container for angle spectrum
    for isc = 1:24   %size( transpose(Scales) )
        anglespec(:,isc) = squeeze(innercoherence(:,isc)) ./ transpose(meanbyscale);
    end 
    
    % Define scales and angles for the upper half
    Angles_pos = linspace(0, pi, size(anglespec, 2));  % Angles from 0 to pi (upper half)
    % Scales = 10.^(1:.05:1.9);  % Logarithmic scales (adjust based on your data)
    
    % Extend angles to cover the full circle (0 to 2*pi)
    Angles_full = [Angles_pos, Angles_pos + pi];  % Full angles from 0 to 2*pi
    
    % Create meshgrid for polar coordinates
    [Theta, R] = meshgrid(Angles_full, Scales);  % Full angle grid for both halves
    
    % Extend the data to cover the lower half, not just mirrored
    anglespec_full = [anglespec, anglespec];  % Repeat the data for the second half (adjust based on physical symmetry)
    
    % Convert polar coordinates to Cartesian for plotting with pcolor
    [X, Y] = pol2cart(Theta, R);
    
    % Plot using Cartesian coordinates with pcolor
    figure;
    pcolor(X, Y, anglespec_full);  % Use pcolor to plot the wavelet power data
    %colormap("sky")
    shading interp;  % Smooth shading
    
    % Customize the plot with labels and title
    title('Wavelet Coherence Full Rose Plot');
    colorbar;  % Add colorbar to indicate wavelet power
    ylabel(colorbar, 'Magnitude');  % Add label to the colorbar
    
    % Remove Cartesian ticks 
    set(gca, 'XTick', [], 'YTick', []);  
    axis equal;  % Ensure equal scaling
    
    % Add angular lines with fractions of pi as labels
    hold on;
    angle_ticks_pos = linspace(0, pi, 7);  % Upper half ticks (0 to pi)
    angle_ticks_neg = linspace(pi, 2*pi, 7); % Lower half ticks (pi to 2*pi)
    angle_labels_pos = {'0', '\pi/6', '\pi/3', '\pi/2', '2\pi/3', '5\pi/6', '\pi'};  % Labels for upper half
    angle_labels_neg = {'\pi', '7\pi/6', '4\pi/3', '3\pi/2', '5\pi/3', '11\pi/6', '2\pi'};  % Labels for lower half
    
    % Display angle ticks for the upper half
    for i = 1:length(angle_ticks_pos)
        angle_rad = angle_ticks_pos(i);
        
        % Upper half
        plot([0 max(X(:))*cos(angle_rad)], [0 max(Y(:))*sin(angle_rad)], 'k--');
        text(max(X(:))*cos(angle_rad)*1.1, max(Y(:))*sin(angle_rad)*1.1, ...
            angle_labels_pos{i}, 'HorizontalAlignment', 'center');
    end
    
    % Display angle ticks for the lower half
    for i = 1:length(angle_ticks_neg)-1
        angle_rad = angle_ticks_neg(i);
        
        % Lower half
        plot([0 max(X(:))*cos(angle_rad)], [0 max(Y(:))*sin(angle_rad)], 'k--');
        text(max(X(:))*cos(angle_rad)*1.1, max(Y(:))*sin(angle_rad)*1.1, ...
            angle_labels_neg{i}, 'HorizontalAlignment', 'center');
    end
    
    % Add radial rings using actual Scales values
    ring_radii = Scales;  % Use actual scales for the rings
    for i = 1:length(ring_radii)
        if ~ismember(i, [2:3 ,5:10,12:13,15,17])
            theta_ring = linspace(0, 2*pi, 100);  % Full circle for each ring
            [x_ring, y_ring] = pol2cart(theta_ring, ring_radii(i));
            plot(x_ring, y_ring, 'k--');  % Plot the radial ring
    
            text(max(x_ring)*1.05, 0, num2str(ring_radii(i), '%.2f'), 'HorizontalAlignment', 'left')
        end
    end
    hold off;
end

function peak_list = plot_waverose_Coherence_Maximums(coherence, Scales)
    buffer = round( max(Scales) );
    % Calculate inner coherence
    innercoherence = squeeze( mean(mean( coherence(buffer:size(coherence,1)-buffer, ...
                                           buffer:size(coherence,2)-buffer, :,:) )));
    meanbyscale = squeeze(mean(transpose(innercoherence)));
    
    % Initialize the angle spectrum container
    anglespec = zeros(size(innercoherence));
    for isc = 1:length(Scales)
        anglespec(:, isc) = squeeze(innercoherence(:, isc)) ./ transpose(meanbyscale);
    end
    
    % Define scales and angles for the upper half
    Angles_pos = linspace(0, pi, size(anglespec, 2));  % Angles from 0 to pi (upper half)
    
    % Extend angles to cover the full circle (0 to 2*pi)
    Angles_full = [Angles_pos, Angles_pos + pi];  % Full angles from 0 to 2*pi
    
    % Create meshgrid for polar coordinates
    [Theta, R] = meshgrid(Angles_full, Scales);  % Full angle grid for both halves
    
    % Extend the data to cover the lower half, not just mirrored
    anglespec_full = [anglespec, anglespec];  % Repeat the data for the second half
    
    % Convert polar coordinates to Cartesian for plotting with pcolor
    [X, Y] = pol2cart(Theta, R);
    
    % Plot using Cartesian coordinates with pcolor
    figure;
    pcolor(X, Y, anglespec_full);  % Use pcolor to plot the wavelet power data
    shading interp;  % Smooth shading
    
    % Customize the plot with labels and title
    title('Wavelet Coherence Full Rose Plot with Local Maxima');
    colorbar;  % Add colorbar to indicate wavelet power
    ylabel(colorbar, 'Magnitude');  % Add label to the colorbar
    
    % Remove Cartesian ticks
    set(gca, 'XTick', [], 'YTick', []);
    axis equal;  % Ensure equal scaling
    
    % Add angular lines with fractions of pi as labels
    hold on;
    angle_ticks_pos = linspace(0, pi, 7);    % Upper half ticks (0 to pi)
    angle_ticks_neg = linspace(pi, 2*pi, 7); % Lower half ticks (pi to 2*pi)
    angle_labels_pos = {'0', '\pi/6', '\pi/3', '\pi/2', '2\pi/3', '5\pi/6', '\pi'};
    angle_labels_neg = {'\pi', '7\pi/6', '4\pi/3', '3\pi/2', '5\pi/3', '11\pi/6', '2\pi'};
    
    % Display angle ticks for the upper half
    for i = 1:length(angle_ticks_pos)
        angle_rad = angle_ticks_pos(i);
        % Upper half
        plot([0, max(X(:))*cos(angle_rad)], [0, max(Y(:))*sin(angle_rad)], 'k--');
        text(max(X(:))*cos(angle_rad)*1.1, max(Y(:))*sin(angle_rad)*1.1, ...
            angle_labels_pos{i}, 'HorizontalAlignment', 'center');
    end
    
    % Display angle ticks for the lower half
    for i = 1:length(angle_ticks_neg)-1
        angle_rad = angle_ticks_neg(i);
        % Lower half
        plot([0, max(X(:))*cos(angle_rad)], [0, max(Y(:))*sin(angle_rad)], 'k--');
        text(max(X(:))*cos(angle_rad)*1.1, max(Y(:))*sin(angle_rad)*1.1, ...
            angle_labels_neg{i}, 'HorizontalAlignment', 'center');
    end
    
    % Add radial rings using actual Scales values
    ring_radii = Scales;  % Use actual scales for the rings
    for i = 1:length(ring_radii)
        if ~ismember(i, [2:3, 5:10, 12:13, 15, 17])
            theta_ring = linspace(0, 2*pi, 100);  % Full circle for each ring
            [x_ring, y_ring] = pol2cart(theta_ring, ring_radii(i));
            plot(x_ring, y_ring, 'k--');  % Plot the radial ring
            
            text(max(x_ring)*1.05, 0, num2str(ring_radii(i), '%.2f'), 'HorizontalAlignment', 'left');
        end
    end
    
    % Step 1: Find local maxima in 'anglespec' (upper half)
    local_maxima = imregionalmax(anglespec);
    
    % Step 2: Exclude edge extrema
    local_maxima(1, :) = false;         % Exclude first row (edge in scale)
    local_maxima(end, :) = false;       % Exclude last row (edge in scale)
    local_maxima(:, 1) = false;         % Exclude first column (edge in angle)
    local_maxima(:, end) = false;       % Exclude last column (edge in angle)
    
    % Step 3: Exclude maxima below a certain threshold
    threshold = 1;  % Set your threshold value here
    local_maxima(anglespec < threshold) = false;
    
    % Step 4: Get indices of remaining maxima
    [row_indices, col_indices] = find(local_maxima);
    
    % Step 5: Map indices to scales and angles
    max_scales = Scales(row_indices);
    max_angles = Angles_pos(col_indices);
    
    % Step 6: Extend maxima to full circle
    max_angles_full = [max_angles; max_angles + pi];
    max_scales_full = [max_scales; max_scales];
    
    % Step 7: Convert maxima to Cartesian coordinates
    [max_X, max_Y] = pol2cart(max_angles_full, max_scales_full);
    
    % Step 8: Overlay maxima on the rose plot
    plot(max_X, max_Y, 'k*', 'MarkerSize', 10);
    
    hold off;

    % Step 9: Prepare the list of scale/angle pairs
    peak_list = [max_scales; max_angles];  % Create a list of scale/angle pairs for the upper half
    
end

function plot_waverose_TrueCoherence(spec1, spec2, Scales, Angles)
    % Compute the individual wavelet power spectra
    power1 = abs(spec1).^2;
    power2 = abs(spec2).^2;
    
    % Compute the cross-wavelet spectrum (XWT)
    xwt = spec1 .* conj(spec2);
    
    % Apply smoothing to the spectra
    % Define the smoothing window size
    window_size = 5;  % Adjust based on your data
    
    % Create a smoothing kernel (2D Gaussian kernel)
    sigma = window_size / 2;
    [x, y] = meshgrid(-window_size:window_size, -window_size:window_size);
    gaussian_kernel = exp(-(x.^2 + y.^2) / (2 * sigma^2));
    gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));  % Normalize the kernel
    
    % Smooth the spectra along spatial dimensions
    smooth_power1 = zeros(size(power1));
    smooth_power2 = zeros(size(power2));
    smooth_xwt = zeros(size(xwt));
    
    % Apply smoothing to each scale and angle
    for s = 1:size(power1, 3)
        for a = 1:size(power1, 4)
            smooth_power1(:,:,s,a) = conv2(power1(:,:,s,a), gaussian_kernel, 'same');
            smooth_power2(:,:,s,a) = conv2(power2(:,:,s,a), gaussian_kernel, 'same');
            smooth_xwt(:,:,s,a) = conv2(xwt(:,:,s,a), gaussian_kernel, 'same');
        end
    end
    
    % Compute the wavelet coherence
    numerator = abs(smooth_xwt).^2;
    denominator = smooth_power1 .* smooth_power2;
    coherence = numerator ./ denominator;
    
    % Cap coherence values at 1
    coherence(coherence > 1) = 1;
    
    % Now proceed to compute the average coherence over a region
    % Adjust the indices as needed
    row_start = 50; row_end = 800;
    col_start = 200; col_end = 800;
    region_coherence = coherence(row_start:row_end, col_start:col_end, :, :);
    
    % Compute the mean coherence over the spatial dimensions (rows and cols)
    innercoherence = squeeze(mean(mean(region_coherence, 1), 2));  % Resulting in [scales x angles]
    
    % Use the coherence values directly
    anglespec = innercoherence;  % [scales x angles]
    
    % Define angles for the upper half
    Angles_pos = Angles;  % Assuming Angles go from 0 to pi
    
    % Extend angles to cover the full circle (0 to 2*pi)
    Angles_full = [Angles_pos, Angles_pos + pi];  % Full angles from 0 to 2*pi
    
    % Extend the data to cover the lower half
    anglespec_full = [anglespec, anglespec];  % Repeat the data for the second half
    
    % Create meshgrid for polar coordinates
    [R, Theta] = meshgrid(Scales, Angles_full);
    
    % Convert polar coordinates to Cartesian for plotting with pcolor
    [X, Y] = pol2cart(Theta', R');  % Transpose to match dimensions
    
    % Plot using Cartesian coordinates with pcolor
    figure;
    pcolor(X, Y, anglespec_full');
    colormap("hot");
    shading interp;  % Smooth shading
    
    % Customize the plot with labels and title
    title('Wavelet Coherence Full Rose Plot');
    colorbar;  % Add colorbar to indicate wavelet coherence
    ylabel(colorbar, 'Wavelet Coherence');  % Add label to the colorbar
    
    % Set color limits between 0 and 1 for coherence
    clim([0, 1]);
    
    % Remove Cartesian ticks
    set(gca, 'XTick', [], 'YTick', []);
    axis equal;  % Ensure equal scaling
    
    % Add angular lines with fractions of pi as labels
    hold on;
    angle_ticks = 0:pi/6:2*pi;  % Ticks every pi/6 from 0 to 2*pi
    angle_labels = {'0', '\pi/6', '\pi/3', '\pi/2', '2\pi/3', '5\pi/6', ...
                    '\pi', '7\pi/6', '4\pi/3', '3\pi/2', '5\pi/3', '11\pi/6', '2\pi'};
    
    % Maximum radius for plotting lines
    max_radius = max(Scales);
    
    % Display angle ticks for the full circle
    for i = 1:length(angle_ticks)
        angle_rad = angle_ticks(i);
        % Plot radial lines
        plot([0 max_radius * cos(angle_rad)], [0 max_radius * sin(angle_rad)], 'k--');
        % Add angle labels
        text(max_radius * cos(angle_rad) * 1.1, max_radius * sin(angle_rad) * 1.1, ...
            angle_labels{i}, 'HorizontalAlignment', 'center');
    end
    
    % Add radial rings using actual Scales values
    ring_radii = Scales;  % Use actual scales for the rings
    for i = 1:length(ring_radii)
        if ~ismember(i, [2:3, 5:10, 12:13, 15, 17])
            theta_ring = linspace(0, 2*pi, 200);  % Full circle for each ring
            [x_ring, y_ring] = pol2cart(theta_ring, ring_radii(i));
            plot(x_ring, y_ring, 'k--');  % Plot the radial ring
            % Place the scale labels
            text(ring_radii(i) * 1.05, 0, num2str(ring_radii(i), '%.2f'), 'HorizontalAlignment', 'left');
        end
    end
    hold off;
end
