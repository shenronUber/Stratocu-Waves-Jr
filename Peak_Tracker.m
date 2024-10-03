% Initialize video reader
v = VideoReader("DATA/closedcellday_2022_09_06.mp4");

% Initialize video writer
outputVideo = VideoWriter('annotated_video.mp4', 'MPEG-4');
open(outputVideo);

% Read total number of frames
num_frames = v.NumFrames; % Total number of frames in the video

% define preprocessing option
Preprocess_Flag = 1; % 1 is on / 0 is off

% Define wavelet transform parameters
NSCALES = 24;
Angles = 0:pi/NSCALES:pi;
Scales = 10.^(1:0.05:1.9);
shrinkfactor = 5; % As per your code
invshrinkfactor = 1 / shrinkfactor;
Scales = Scales / shrinkfactor;

% Time interval between frames (adjust if necessary)
time_interval = 1800; % Assuming 30 minutes in seconds

% Pixel size in kilometers (adjust based on your data)
pixel_size_km = 2 * shrinkfactor;

% Initialize variables to store peak locations over time
% We'll collect all peak information in a struct array
all_peaks = [];

% Loop over frame pairs
for frame_idx = 1:(num_frames - 1)
    % Read frames
    frame1 = read(v, frame_idx);
    frame2 = read(v, frame_idx + 1);

    % Extract the red channel and convert to double precision
    red1 = double(frame1(:, :, 1));
    red2 = double(frame2(:, :, 1));
    
    %preprocess frames
    if Preprocess_Flag == 1 
        red1 = preprocess_img(red1);
        red2 = preprocess_img(red2);
    end

    % Resize images to make computations faster
    frame1_resized = imresize(red1, invshrinkfactor);
    frame2_resized = imresize(red2, invshrinkfactor);

    % Get dimensions of resized frames
    [h, w] = size(frame1_resized);

    % Perform 2D Continuous Wavelet Transform on both frames
    cwt1 = cwtft2(frame1_resized, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
    spec1 = squeeze(cwt1.cfs);

    cwt2 = cwtft2(frame2_resized, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
    spec2 = squeeze(cwt2.cfs);

    % Compute the cross-wavelet spectrum (XWT)
    xwt = spec1 .* conj(spec2);

    % Compute coherence and phase difference
    power1 = abs(spec1).^2;
    power2 = abs(spec2).^2;
    coherence = abs(xwt).^2 ./ (power1 .* power2);
    phase_difference = angle(xwt);

    % Compute area-averaged power and normalize by mean increase with scale
    buffer = round(max(Scales));
    innerpower = squeeze(mean(mean(power1(buffer:end-buffer, buffer:end-buffer, :, :))));
    meanbyscale = squeeze(mean(transpose(innerpower)));

    % Normalize by the mean increase with scale
    anglespec = innerpower .* 0;
    for isc = 1:length(Scales)
        anglespec(:, isc) = squeeze(innerpower(:, isc)) ./ transpose(meanbyscale);
    end

    % Find peaks where anglespec exceeds a threshold strength
    peak_threshold = 1.8; % Increased threshold to 80%
    [row, col] = find(imregionalmax(anglespec) & anglespec > peak_threshold);

    % Initialize a struct array to store peaks for this frame
    frame_peaks = [];

    % For each peak in anglespec
    for peak_idx = 1:length(row)
        isc = row(peak_idx);
        ian = col(peak_idx);

        % Get the scale and angle corresponding to the peak
        scale = Scales(isc);
        ang = Angles(ian);

        % Extract the corresponding power map for the peak
        power_map = abs(spec1(:, :, isc, ian)).^2;

        % Find spatial peaks in the power map
        spatial_peaks = imregionalmax(power_map);
        power_threshold_map = 0.8 * max(power_map(:)); % Increased to 80%
        spatial_peaks = spatial_peaks & (power_map >= power_threshold_map);

        % Get coordinates of spatial peaks
        [y_peaks, x_peaks] = find(spatial_peaks);

        % Ensure there are peaks detected
        if isempty(x_peaks)
            continue; % Skip to the next peak if none are found
        end

        % Calculate the mean phase difference over spatial peaks
        phase_slice = phase_difference(:, :, isc, ian);
        phase_mean = mean(phase_slice(spatial_peaks));

        % Calculate speed
        wavelength_km = scale * 2 * pixel_size_km; % Since scale is half the wavelength
        distance_shift_km = phase_mean * wavelength_km / (2 * pi);
        speed_km_per_s = distance_shift_km / time_interval;

        % Calculate propagation direction (angle + 90 degrees)
        propagation_direction = ang + pi / 2;

        % Map coordinates back to original frame size
        x_peaks_full = x_peaks * shrinkfactor;
        y_peaks_full = y_peaks * shrinkfactor;

        % Store peak information in struct array
        for i = 1:length(x_peaks_full)
            peak_info.frame_idx = frame_idx;
            peak_info.x = x_peaks_full(i);
            peak_info.y = y_peaks_full(i);
            peak_info.scale = scale;
            peak_info.angle = ang;
            peak_info.phase_mean = phase_mean;
            peak_info.speed = speed_km_per_s * 1000; % Convert to m/s
            peak_info.propagation_direction = propagation_direction;
            frame_peaks = [frame_peaks; peak_info];
        end
    end

    % Append peaks from this frame to all_peaks
    all_peaks = [all_peaks; frame_peaks];
end

%% Tracking Algorithm Based on Position and Scale/Angle Index Matching
%{
% Define tolerances in index steps
scale_tolerance_steps = 1; % Allow scales within +/- 1 index
angle_tolerance_steps = 1; % Allow angles within +/- 1 index

% Get total number of frames
num_frames = max([all_peaks.frame_idx]);

% Initialize peaks_by_frame as an empty cell array
peaks_by_frame = cell(num_frames, 1);

% Populate peaks_by_frame
for i = 1:length(all_peaks)
    frame_idx = all_peaks(i).frame_idx;
    peaks_by_frame{frame_idx} = [peaks_by_frame{frame_idx}, i]; % Store index of peak in that frame
end

% Load the first frame and calculate its size
frame1 = read(v, 1); % Read the first frame of the video
image_size = size(frame1(:,:,1)); % Assume it's the same size for all frames
max_distance = 0.25 * sqrt(image_size(1)^2 + image_size(2)^2)*0.25; % Max movement is a quarter of diagonal

% Initialize tracking parameters
next_track_id = 1;

% Assign track_id = 0 to all peaks initially, and compute scale and angle indices
for i = 1:length(all_peaks)
    all_peaks(i).track_id = 0;
    % Find indices in Scales and Angles arrays
    all_peaks(i).scale_idx = find(abs(Scales - all_peaks(i).scale) < 1e-6, 1);
    all_peaks(i).angle_idx = find(abs(Angles - all_peaks(i).angle) < 1e-6, 1);
end

% For frame 1, assign new track IDs to all peaks
for p_idx = peaks_by_frame{1}
    all_peaks(p_idx).track_id = next_track_id;
    next_track_id = next_track_id + 1;
end

% Loop over frame pairs to link peaks based on position and scale/angle index matching
for n = 1:(num_frames - 1)
    peaks_n = peaks_by_frame{n};
    peaks_n1 = peaks_by_frame{n + 1};
    
    % For each peak in frame n
    for p_idx = peaks_n
        track_id = all_peaks(p_idx).track_id;
        x = all_peaks(p_idx).x;
        y = all_peaks(p_idx).y;
        scale_idx = all_peaks(p_idx).scale_idx;
        angle_idx = all_peaks(p_idx).angle_idx;

        % Find peaks in frame n+1 that are within max_distance, not yet assigned, and share scale or angle within tolerance
        candidate_indices = [];
        distances = [];
        for q_idx = peaks_n1
            if all_peaks(q_idx).track_id == 0
                x_q = all_peaks(q_idx).x;
                y_q = all_peaks(q_idx).y;
                scale_idx_q = all_peaks(q_idx).scale_idx;
                angle_idx_q = all_peaks(q_idx).angle_idx;
                
                % Check if scales or angles match within tolerance steps
                scales_match = abs(scale_idx_q - scale_idx) <= scale_tolerance_steps;
                angles_match = abs(angle_idx_q - angle_idx) <= angle_tolerance_steps;
                
                if scales_match || angles_match
                    % Calculate the distance between the two points
                    dx = x_q - x;
                    dy = y_q - y;
                    distance = sqrt(dx^2 + dy^2);

                    % Ensure the distance is below the maximum allowable distance
                    if distance <= max_distance
                        candidate_indices = [candidate_indices, q_idx];
                        distances = [distances, distance];
                    end
                end
            end
        end
        
        if ~isempty(candidate_indices)
            % Find the closest candidate based on spatial distance
            [~, min_idx] = min(distances);
            matched_q_idx = candidate_indices(min_idx);
            
            % Assign the same track_id to the matched peak
            all_peaks(matched_q_idx).track_id = track_id;
        end
    end

    % Now assign new track IDs to any unmatched peaks in frame n+1
    for q_idx = peaks_n1
        if all_peaks(q_idx).track_id == 0
            % Assign a new track ID
            all_peaks(q_idx).track_id = next_track_id;
            next_track_id = next_track_id + 1;
        end
    end
end

% Collect tracks by grouping peaks with the same track_id
unique_track_ids = unique([all_peaks.track_id]);
unique_track_ids(unique_track_ids == 0) = []; % Remove unassigned peaks (track_id == 0)

tracks = struct('track_id', {}, 'frames', {}, 'x', {}, 'y', {}, 'scales', {}, 'angles', {});

for i = 1:length(unique_track_ids)
    track_id = unique_track_ids(i);
    % Get indices of peaks with this track_id
    track_peak_indices = find([all_peaks.track_id] == track_id);
    % Sort by frame index
    [~, sort_idx] = sort([all_peaks(track_peak_indices).frame_idx]);
    track_peak_indices = track_peak_indices(sort_idx);
    % Get frames
    frames = [all_peaks(track_peak_indices).frame_idx];
    % Get x, y
    x = [all_peaks(track_peak_indices).x];
    y = [all_peaks(track_peak_indices).y];
    % Get scales and angles (may vary over time now)
    scales = [all_peaks(track_peak_indices).scale];
    angles = [all_peaks(track_peak_indices).angle];
    % Store in tracks
    tracks(i).track_id = track_id;
    tracks(i).frames = frames;
    tracks(i).x = x;
    tracks(i).y = y;
    tracks(i).scales = scales;
    tracks(i).angles = angles;
end

% Filter tracks that span at least 3 frames
valid_tracks = tracks(arrayfun(@(t) length(t.frames) >= 3, tracks));

%}

%% Implement Improved Tracking Algorithm Based on Cost Function (Position, Scale, and Angle)

% Define the weights for the cost function
w_position = 1; % Weight for position difference
w_scale = 0.5;  % Weight for scale difference
w_angle = 0.5;  % Weight for angle difference

% Define the maximum allowable cost
max_cost = 1.5; % Adjust this threshold based on empirical observation

% Get total number of frames
num_frames = max([all_peaks.frame_idx]);

% Initialize peaks_by_frame as an empty cell array
peaks_by_frame = cell(num_frames, 1);

% Populate peaks_by_frame
for i = 1:length(all_peaks)
    frame_idx = all_peaks(i).frame_idx;
    peaks_by_frame{frame_idx} = [peaks_by_frame{frame_idx}, i]; % Store index of peak in that frame
end

% Load the first frame and calculate its size
frame1 = read(v, 1); % Read the first frame of the video
image_size = size(frame1(:,:,1)); % Assume it's the same size for all frames
max_distance = (0.25/2) * sqrt(image_size(1)^2 + image_size(2)^2); % Max movement is half a quarter of diagonal

% Initialize tracking parameters
next_track_id = 1;

% Assign track_id = 0 to all peaks initially
for i = 1:length(all_peaks)
    all_peaks(i).track_id = 0;
end

% For frame 1, assign new track IDs to all peaks
for p_idx = peaks_by_frame{1}
    all_peaks(p_idx).track_id = next_track_id;
    next_track_id = next_track_id + 1;
end

% Cost function calculation
calculate_cost = @(dx, dy, scale_diff, angle_diff, max_distance) ...
    (w_position * sqrt(dx^2 + dy^2) / max_distance) + ...
    (w_scale * abs(scale_diff)) + ...
    (w_angle * abs(angle_diff));

% Loop over frame pairs to link peaks based on the cost function
for n = 1:(num_frames - 1)
    peaks_n = peaks_by_frame{n};
    peaks_n1 = peaks_by_frame{n + 1};
    
    % For each peak in frame n
    for p_idx = peaks_n
        track_id = all_peaks(p_idx).track_id;
        x = all_peaks(p_idx).x;
        y = all_peaks(p_idx).y;
        scale = all_peaks(p_idx).scale;
        angle = all_peaks(p_idx).angle;

        % Initialize arrays to store candidate indices and costs
        candidate_indices = [];
        candidate_costs = [];

        % For each peak in frame n+1
        for q_idx = peaks_n1
            if all_peaks(q_idx).track_id == 0
                x_q = all_peaks(q_idx).x;
                y_q = all_peaks(q_idx).y;
                scale_q = all_peaks(q_idx).scale;
                angle_q = all_peaks(q_idx).angle;

                % Compute differences
                dx = x_q - x;
                dy = y_q - y;
                scale_diff = scale_q - scale;
                angle_diff = min(abs(angle_q - angle), 2*pi - abs(angle_q - angle)); % Consider periodicity of angle

                % Compute the total cost
                cost = calculate_cost(dx, dy, scale_diff, angle_diff, max_distance);

                % If the cost is below the maximum allowable cost, consider it a candidate
                if cost <= max_cost
                    candidate_indices = [candidate_indices, q_idx];
                    candidate_costs = [candidate_costs, cost];
                end
            end
        end
        
        if ~isempty(candidate_indices)
            % Find the candidate with the minimum cost
            [~, min_idx] = min(candidate_costs);
            matched_q_idx = candidate_indices(min_idx);
            
            % Assign the same track_id to the matched peak
            all_peaks(matched_q_idx).track_id = track_id;
        end
    end

    % Now assign new track IDs to any unmatched peaks in frame n+1
    for q_idx = peaks_n1
        if all_peaks(q_idx).track_id == 0
            % Assign a new track ID
            all_peaks(q_idx).track_id = next_track_id;
            next_track_id = next_track_id + 1;
        end
    end
end

% Collect tracks by grouping peaks with the same track_id
unique_track_ids = unique([all_peaks.track_id]);
unique_track_ids(unique_track_ids == 0) = []; % Remove unassigned peaks (track_id == 0)

tracks = struct('track_id', {}, 'frames', {}, 'x', {}, 'y', {}, 'scales', {}, 'angles', {});

for i = 1:length(unique_track_ids)
    track_id = unique_track_ids(i);
    % Get indices of peaks with this track_id
    track_peak_indices = find([all_peaks.track_id] == track_id);
    % Sort by frame index
    [~, sort_idx] = sort([all_peaks(track_peak_indices).frame_idx]);
    track_peak_indices = track_peak_indices(sort_idx);
    % Get frames
    frames = [all_peaks(track_peak_indices).frame_idx];
    % Get x, y
    x = [all_peaks(track_peak_indices).x];
    y = [all_peaks(track_peak_indices).y];
    % Get scales and angles (may vary over time now)
    scales = [all_peaks(track_peak_indices).scale];
    angles = [all_peaks(track_peak_indices).angle];
    % Store in tracks
    tracks(i).track_id = track_id;
    tracks(i).frames = frames;
    tracks(i).x = x;
    tracks(i).y = y;
    tracks(i).scales = scales;
    tracks(i).angles = angles;
end

% Filter tracks that span at least 3 frames
valid_tracks = tracks(arrayfun(@(t) length(t.frames) >= 3, tracks));


%% Visualize Tracks Over Time with Updated Labels

% Initialize video writer
output_video = VideoWriter('tracked_peaks_video.avi'); % Name of the output video file
output_video.FrameRate = 10; % Frames per second for the video
open(output_video); % Open the video writer

% Assign a unique color to each track
num_tracks = length(valid_tracks);
colors = lines(num_tracks);

% Loop through the frames and overlay the tracked peaks
for frame_idx = 1:num_frames
    % Read the current frame
    current_frame = read(v, frame_idx);
    
    % Display the frame
    figure(1); clf;
    imshow(current_frame);
    hold on;

    % For each valid track, plot the trajectory up to the current frame
    for i = 1:num_tracks
        track_frames = valid_tracks(i).frames;
        % Check if the current frame is part of the track
        if ismember(frame_idx, track_frames)
            % Find the index of the current frame in the track
            frame_in_track_idx = find(track_frames == frame_idx);
            % Get the x, y position of the peak in this frame
            x = valid_tracks(i).x(frame_in_track_idx);
            y = valid_tracks(i).y(frame_in_track_idx);
            % Get the scale and angle at this frame
            scale = valid_tracks(i).scales(frame_in_track_idx);
            angle = valid_tracks(i).angles(frame_in_track_idx);
            % Plot the peak's current position
            plot(x, y, 'o', 'Color', colors(i, :), 'MarkerSize', 5, 'LineWidth', 2);
            % Annotate with scale and angle at this frame
            text(x, y, sprintf('Scale: %.2f\nAngle: %.1f°', scale, rad2deg(angle)), ...
                 'Color', colors(i, :), 'FontSize', 10);
            
            % Plot the entire trajectory of the track up to this point
            plot(valid_tracks(i).x(1:frame_in_track_idx), valid_tracks(i).y(1:frame_in_track_idx), ...
                 '-', 'Color', colors(i, :), 'LineWidth', 2);
        end
    end

    title(sprintf('Tracked Peaks - Frame %d', frame_idx));
    hold off;
    
    % Capture the frame for the video
    frame = getframe(gcf); 
    writeVideo(output_video, frame); % Write the frame to the video
end

% Close the video writer
close(output_video);

disp('Video creation complete! The tracked peaks video has been saved.');

%% Static Visualization of Tracks on the First Frame

% Read the first frame for background
frame1 = read(v, 1);
figure;
imshow(frame1);
hold on;

% For each valid track, plot the entire trajectory
for i = 1:num_tracks
    % Plot the trajectory
    plot(valid_tracks(i).x, valid_tracks(i).y, '-o', 'Color', colors(i, :), 'LineWidth', 2, 'MarkerSize', 5);
    % Annotate with initial scale and angle
    text(valid_tracks(i).x(1), valid_tracks(i).y(1), sprintf('Scale: %.2f\nAngle: %.1f°', ...
        valid_tracks(i).scales(1), rad2deg(valid_tracks(i).angles(1))), 'Color', colors(i, :), 'FontSize', 10);
end

title('Tracked Peaks Over Time');
hold off;

%% Visualise single frames 

% User selects which frame to plot
frame_to_plot = 1;  % Adjust this value to choose a specific frame

% Load the frame image
frame_image = read(v, frame_to_plot);
if Preprocess_Flag == 1
    frame_image_resized = imresize(double(preprocess_img(frame_image(:,:,1))), invshrinkfactor);  % Adjust for the resize
else
    frame_image_resized = imresize(double(frame_image(:,:,1)), invshrinkfactor);  % Adjust for the resize
end
% Extract the wavelet spectrum for the chosen frame
cwt = cwtft2(frame_image_resized, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
spec = squeeze(cwt.cfs);

% Calculate power spectrum
power_spec = abs(spec).^2;

% Extract peaks for the selected frame
selected_peak_indices = peaks_by_frame{frame_to_plot};  % Get indices of peaks
selected_peaks = all_peaks(selected_peak_indices);      % Retrieve actual peak structs

% Initialize figure for subplots
figure;
tiledlayout(1, 2);

% Assign colors based on unique combinations of scale and angle
unique_peak_params = unique([[selected_peaks.scale]', [selected_peaks.angle]'], 'rows');
num_unique_peaks = size(unique_peak_params, 1);
colors = lines(num_unique_peaks);  % Assign distinct colors for each unique scale-angle pair

% Create a map to store the color for each scale-angle pair
peak_color_map = containers.Map('KeyType', 'char', 'ValueType', 'any');

% Assign colors to each unique scale-angle pair
for i = 1:num_unique_peaks
    scale_angle_str = sprintf('%.6f_%.6f', unique_peak_params(i, 1), unique_peak_params(i, 2));
    peak_color_map(scale_angle_str) = colors(i, :);  % Store the color for this unique scale-angle pair
end

% Plot 1: Power Rose Spectrum with Peaks and Speed
nexttile;

% Compute the mean inner power within the buffered region
inner_power = squeeze(mean(mean(power_spec(buffer:size(power_spec, 1)-buffer, ...
                                          buffer:size(power_spec, 2)-buffer, :, :))));

mean_by_scale = squeeze(mean(transpose(inner_power)));

% Normalize by the mean increase with scale
angle_spec = inner_power .* 0;
for isc = 1:length(Scales)
    angle_spec(:, isc) = squeeze(inner_power(:, isc)) ./ transpose(mean_by_scale);
end

% Define angles for the upper half
Angles_pos = Angles;  % Use the Angles vector defined earlier

% Ensure that X, Y, and angle_spec have consistent dimensions
[Theta, R] = meshgrid(Angles_pos, Scales);  % Matching dimensions of Scales and Angles

% Convert polar coordinates to Cartesian for plotting
[X, Y] = pol2cart(Theta, R);

% Plot using Cartesian coordinates with pcolor
pcolor(X, Y, angle_spec);
shading interp;
colormap('jet');
colorbar;
title('Power Rose Spectrum with Peaks and Speed');
axis equal
axis([-16 16 0 20])
hold on;


% Overlay detected peaks with their speed, using the synchronized colors
fprintf('Peaks in Rose Plot:\n');
for i = 1:length(selected_peaks)
    % Get the scale and angle of the current peak
    scale = selected_peaks(i).scale;
    angle = selected_peaks(i).angle;

    % Use the scale and angle to find the corresponding color
    scale_angle_str = sprintf('%.6f_%.6f', scale, angle);
    peak_color = peak_color_map(scale_angle_str);

    % Plot the peak with the synchronized color
    [x_peak, y_peak] = pol2cart(angle, scale);
    plot(x_peak, y_peak, 'o', 'MarkerSize', 10, 'MarkerFaceColor', peak_color, 'MarkerEdgeColor', 'k');

    % Annotate with speed
    text(x_peak, y_peak, sprintf('%.2f m/s', selected_peaks(i).speed), ...
        'Color', peak_color, 'FontSize', 10, 'HorizontalAlignment', 'left');

    % Print peak details for debugging
    fprintf('Peak %d: Scale = %.2f, Angle = %.1f°, Speed = %.2f m/s\n', ...
            i, scale, rad2deg(angle), selected_peaks(i).speed);
end

hold off;

% Plot 2: Real Image with Overlayed Positions of Power Spectrum
nexttile;

imshow(frame_image, []);
hold on;

% Overlay detected peak positions with angle and scale, using synchronized colors
fprintf('Peaks in Real Image:\n');
for i = 1:length(selected_peaks)
    % Get the position, scale, and angle of the peak
    x = selected_peaks(i).x;
    y = selected_peaks(i).y;
    scale = selected_peaks(i).scale;
    angle = selected_peaks(i).angle;

    % Use the scale and angle to find the corresponding color
    scale_angle_str = sprintf('%.6f_%.6f', scale, angle);
    peak_color = peak_color_map(scale_angle_str);

    % Plot with the synchronized color for each blob
    plot(x, y, 'o', 'MarkerSize', 10, 'MarkerFaceColor', peak_color, 'MarkerEdgeColor', 'k');

    % Annotate with angle and scale
    text(x, y, sprintf('Angle: %.1f°\nScale: %.2f', rad2deg(angle), scale), ...
        'Color', peak_color, 'FontSize', 10, 'HorizontalAlignment', 'left');

    % Print peak details for debugging
    fprintf('Peak %d: x = %.2f, y = %.2f, Scale = %.2f, Angle = %.1f°\n', ...
            i, x, y, scale, rad2deg(angle));
end

title('Real Image with Overlayed Power Spectrum Positions');
hold off;

% We'll generate a new figure for each unique scale-angle combination
for i = 1:num_unique_peaks
    % Get the current unique scale and angle
    current_scale = unique_peak_params(i, 1);
    current_angle = unique_peak_params(i, 2);

    % Find the scale and angle indices in Scales and Angles
    [~, scale_idx] = min(abs(Scales - current_scale));
    [~, angle_idx] = min(abs(Angles - current_angle));
    
    % Set contour level factor for better visualization
    clevfactor = 2; % Adjust this value as needed for visibility

    % Create a new figure for this specific (scale, angle) pair
    figure;
    image_with_wavelet_overlay(frame_image_resized, spec, Scales, scale_idx, angle_idx, clevfactor,Preprocess_Flag);
    hold on;

    % Plot all points corresponding to this scale and angle on the figure
    for j = 1:length(selected_peaks)
        if abs(selected_peaks(j).scale - current_scale) < 1e-6 && abs(selected_peaks(j).angle - current_angle) < 1e-6
            % Plot the point on the wavelet overlay
            plot(selected_peaks(j).x * invshrinkfactor, selected_peaks(j).y * invshrinkfactor, 'o', 'MarkerSize', 10, ...
                 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', 'k');
            
            % Annotate with angle and scale for clarity
            text(selected_peaks(j).x * invshrinkfactor, selected_peaks(j).y * invshrinkfactor, ...
                sprintf('x=%.1f, y=%.1f', selected_peaks(j).x * invshrinkfactor, selected_peaks(j).y * invshrinkfactor), ...
                'Color', colors(i, :), 'FontSize', 10, 'HorizontalAlignment', 'left');
        end
    end

    % Add title and close the hold
    title(sprintf('Overlay for Scale: %.2f, Angle: %.1f°', current_scale, rad2deg(current_angle)));
    hold off;
end


%%

function image_with_wavelet_overlay(img, spec, Scales, scale_idx, angle_idx, clevfactor,ProcessFlag)
    if ProcessFlag ==1
        % Normalize img for display
        img_display = img - min(img(:));  % Shift so that minimum is zero
        img_display = img_display / max(img_display(:));  % Scale to [0,1]

        imagesc(img_display); colormap(gray); colorbar; axis on
        hold on

        % Extract the wavelet coefficients at the specified scale and angle
        wavelet_real = real(spec(:, :, scale_idx, angle_idx));
        wavelet_abs = abs(spec(:, :, scale_idx, angle_idx));

        % Adjust contour levels based on the data range
        clevfactor = clevfactor/1.9;
        
        max_real = max(abs(wavelet_real(:))) / clevfactor;
        max_abs = max(wavelet_abs(:)) / clevfactor;

        % Set contour levels for real part
        posLevels = linspace(0.1 * max_real, max_real, 5);
        negLevels = -posLevels;

        % Plot contours of the positive real part
        contour(wavelet_real, 'LevelList', posLevels, 'LineColor', 'red', 'LineWidth', 1);

        % Plot contours of the negative real part
        contour(wavelet_real, 'LevelList', negLevels, 'LineColor', 'blue', 'LineWidth', 1);

        % Plot contours of the power (magnitude squared)
        power_levels = linspace(0.1 * max_abs^2, max_abs^2, 5);
        contour(wavelet_abs.^2, 'LevelList', power_levels, 'LineColor', 'white', 'LineWidth', 1);

        hold off;
    else
        % Overlay wavelet power on image
        image(img); colormap(gray); colorbar; axis on
        hold on

        posLevels = (1:2:9) / clevfactor;
        negLevels = (-9:2:-1) / clevfactor;

        % Adjust contour levels by scale as factor
        sfactor = Scales(scale_idx);

        % Real part is crests and troughs
        contour(real(spec(:, :, scale_idx, angle_idx)), 'LevelList', posLevels * sfactor, 'EdgeColor', 'red');
        contour(real(spec(:, :, scale_idx, angle_idx)), 'LevelList', negLevels * sfactor, 'EdgeColor', 'blue');

        % Power is magnitude squared
        contour((abs(spec(:, :, scale_idx, angle_idx)) * sfactor) .^2, 'EdgeColor', 'white');
        hold off;
    end
end
