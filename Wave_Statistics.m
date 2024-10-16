%% Parameters and Setup

% Start timing the entire process
total_time = tic; % Start the overall timer

% Video file
video_file = "DATA/closedcellday_2022_09_06.mp4";

% Shrinking factor
shrinkfactor = 5;
invshrinkfactor = 1 / shrinkfactor;

% Dynamic pixel size (km per pixel)
original_pixel_size_km = 2; % Original pixel size before shrinking
pixel_size_km = original_pixel_size_km * shrinkfactor; % Adjusted pixel size due to shrinking

% Square size in degrees
square_size_deg = 5; % 5x5 degrees squares

% Conversion factor: 1 degree ≈ 111.32 km on Earth's surface
km_per_degree = 111.32;

% Calculate square size in km
square_size_km = square_size_deg * km_per_degree; % Total km per square

% Calculate square size in pixels
square_size_px = round(square_size_km / pixel_size_km);

% Brightness thresholds
brightness_threshold = 0; % Mean brightness below which squares are ignored 0.4;
std_threshold = 10;      % Standard deviation above which squares are ignored 0.16; 

% Wavelet transform parameters
NSCALES = 24;
Angles = 0:pi/NSCALES:pi;
% Scales ranging from 10 km to 100 km
min_scale_km = 10;
max_scale_km = 100;  % Updated as per your new value
Scales_km = logspace(log10(min_scale_km), log10(max_scale_km), 10);
% Adjust scales to pixels
Scales = Scales_km / pixel_size_km;

% Windowing parameters
window_buffer = round(10 / shrinkfactor); % Adjust window buffer based on shrink factor
window_function = @hann; % Window function to apply (Hann window)

% Preprocessing flag
Preprocess_Flag = 1; % 1 is on / 0 is off

% Time interval between frames (adjust if necessary)
time_interval = 1800; % Assuming 30 minutes in seconds

%% Read Video and Initialize

% Initialize video reader
v = VideoReader(video_file);

% Read total number of frames
num_frames = v.NumFrames; % Total number of frames in the video

% Read the first frame to get dimensions
frame1 = read(v, 1);

% Resize the frame according to the shrink factor
frame1 = imresize(frame1, invshrinkfactor);

[frame_height, frame_width, ~] = size(frame1);

% Calculate number of squares along width and height
num_squares_x = ceil(frame_width / square_size_px);
num_squares_y = ceil(frame_height / square_size_px);

% Generate Squares Coordinates

% Generate coordinates for squares
squares = [];
idx = 1;
for i = 1:num_squares_y
    for j = 1:num_squares_x
        % Calculate pixel indices
        x_start = floor((j - 1) * frame_width / num_squares_x) + 1;
        y_start = floor((i - 1) * frame_height / num_squares_y) + 1;
        x_end = floor(j * frame_width / num_squares_x);
        y_end = floor(i * frame_height / num_squares_y);

        % Apply window buffer
        x_start = x_start + window_buffer;
        y_start = y_start + window_buffer;
        x_end = x_end - window_buffer;
        y_end = y_end - window_buffer;

        % Ensure valid ranges
        if x_end > x_start && y_end > y_start
            squares(idx).x_range = x_start:x_end;
            squares(idx).y_range = y_start:y_end;
            squares(idx).index = idx;
            idx = idx + 1;
        end
    end
end

% Display Frame 1 with Squares

figure;
imshow(frame1);
hold on;
for idx = 1:length(squares)
    rectangle('Position', [squares(idx).x_range(1), squares(idx).y_range(1), ...
        length(squares(idx).x_range), length(squares(idx).y_range)], ...
        'EdgeColor', 'r', 'LineWidth', 1);
    % Label the square
    text(squares(idx).x_range(1), squares(idx).y_range(1) - 5, ...
        sprintf('%d', squares(idx).index), 'Color', 'yellow', 'FontSize', 8);
end
title('Frame 1 with Squares');
hold off;

%% Process Frames and Squares

% Estimate maximum number of peaks per square (adjust based on expected data)
max_peaks_per_square = 10;

% Calculate total number of entries (assuming every square in every frame has max peaks)
total_entries = (num_frames - 1) * length(squares) * max_peaks_per_square;

% Initialize a big NaN array to store results
% Columns: frame_idx, square_idx, scale, angle, phase_difference, speed
peak_data = NaN(total_entries, 6);

% Initialize a counter for the number of entries
entry_counter = 0;

% Measure the time taken for the frame processing loop
frame_processing_time = tic; % Start timing the frame processing loop

% Loop over frames
for frame_idx = 1:num_frames - 1
    fprintf('Processing frame %d/%d...\n', frame_idx, num_frames);

    % Measure time for processing each frame (optional)
    frame_time = tic; % Start timing for a specific frame

    % Read frames
    frame1 = read(v, frame_idx);
    frame2 = read(v, frame_idx + 1);

    % Resize frames according to shrink factor
    frame1 = imresize(frame1, invshrinkfactor);
    frame2 = imresize(frame2, invshrinkfactor);

    % Convert to grayscale if necessary and extract red channel
    if size(frame1, 3) > 1
        frame1 = double(frame1(:, :, 1));
        frame2 = double(frame2(:, :, 1));
    else
        frame1 = double(frame1);
        frame2 = double(frame2);
    end

    % Preprocess frames if needed
    if Preprocess_Flag
        frame1 = preprocess_img(frame1);
        frame2 = preprocess_img(frame2);
    end

    % Loop over squares
    for s_idx = 1:length(squares)
        square = squares(s_idx);

        % Extract the square from both frames
        img1 = frame1(square.y_range, square.x_range);
        img2 = frame2(square.y_range, square.x_range);

        % Check brightness and standard deviation
        mean_brightness = mean(img1(:));
        std_brightness = std(img1(:));

        if mean_brightness < brightness_threshold || std_brightness > std_threshold
            % Ignore this square
            continue;
        end

        % Apply windowing to the square to reduce edge effects
        img1_windowed = apply_window(img1, window_function);
        img2_windowed = apply_window(img2, window_function);

        % Perform wavelet analysis
        img1_windowed = double(img1_windowed);
        img2_windowed = double(img2_windowed);

        % Perform 2D Continuous Wavelet Transform on both frames
        cwt1 = cwtft2(img1_windowed, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
        spec1 = squeeze(cwt1.cfs);

        cwt2 = cwtft2(img2_windowed, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
        spec2 = squeeze(cwt2.cfs);

        % Compute the cross-wavelet spectrum (XWT)
        xwt = spec1 .* conj(spec2);

        % Compute coherence and phase difference
        power1 = abs(spec1).^2;
        power2 = abs(spec2).^2;
        coherence = abs(xwt);
        phase_difference = angle(xwt);
 

        % Find peaks in power spectrum and calculate speeds
        peak_list = find_peaks_and_speeds(coherence, phase_difference, Scales, Angles, pixel_size_km, time_interval);

        % Store the data in the numerical array
        n_peaks = size(peak_list, 2); % Number of peaks found

        if n_peaks > 0
            % Prepare data to store
            data_to_store = [repmat(frame_idx, n_peaks, 1), repmat(s_idx, n_peaks, 1), ...
                             peak_list(1, :)', peak_list(2, :)', peak_list(3, :)', peak_list(4, :)'];

            % Update entry counter
            idx_start = entry_counter + 1;
            idx_end = entry_counter + n_peaks;

            % Ensure we do not exceed preallocated array size
            if idx_end > total_entries
                % Expand the array if necessary
                extra_entries = total_entries; % Double the size
                peak_data = [peak_data; NaN(extra_entries, 6)];
                total_entries = total_entries + extra_entries;
            end

            % Store data
            peak_data(idx_start:idx_end, :) = data_to_store;

            % Update the entry counter
            entry_counter = idx_end;
        end
    end
    
    % Display time taken to process this frame
    fprintf('Time to process frame %d: %.2f seconds.\n', frame_idx, toc(frame_time));
end

% Display total time taken to process all frames
fprintf('Total frame processing time: %.2f seconds.\n', toc(frame_processing_time));

% Remove unused (NaN) entries from peak_data
peak_data = peak_data(1:entry_counter, :);

% Display the overall time taken for the entire routine
fprintf('Total execution time: %.2f seconds.\n', toc(total_time));

%% Analyse the data

analyze_peak_data(peak_data, num_squares_x, num_squares_y)
%plot_waverose(1,20)

%% Supporting Functions

function img_processed = preprocess_img(img)
    % Preprocess the image (e.g., normalization)
    img_processed = img;
    % Normalize to [0, 1]
    img_processed = (img_processed - min(img_processed(:))) / (max(img_processed(:)) - min(img_processed(:)));
end

function img_windowed = apply_window(img, window_function)
    % Apply windowing to reduce edge effects
    [rows, cols] = size(img);
    win_row = window_function(rows);
    win_col = window_function(cols);
    window = win_row * win_col';
    img_windowed = img .* window;
end

function peak_list = find_peaks_and_speeds(coherence, phase_difference, Scales, Angles, pixel_size_km, time_interval)

    buffer = round( max(Scales) );
    % Calculate inner coherence
    innercoherence = squeeze( mean(mean( coherence(buffer:size(coherence,1)-buffer, ...
                                           buffer:size(coherence,2)-buffer, :,:) )));
    meanbyscale = squeeze(mean(transpose(innercoherence)));
    
    % Initialize the angle spectrum container
    anglespec = zeros(size(innercoherence));
    for isc = 1:length(Angles)
        anglespec(:, isc) = squeeze(innercoherence(:, isc)) ./ transpose(meanbyscale);
    end
  
    Angles_pos = linspace(0, pi, size(anglespec, 2));

    % Step 1: Find local maxima in 'anglespec' (upper half)
    local_maxima = imregionalmax(anglespec);
    
    % Step 2: Exclude edge extrema
    %local_maxima(1, :) = false;         % Exclude first row (edge in scale)
    %local_maxima(end, :) = false;       % Exclude last row (edge in scale)
    %local_maxima(:, 1) = false;         % Exclude first column (edge in angle)
    %local_maxima(:, end) = false;       % Exclude last column (edge in angle)
    
    % Step 3: Exclude maxima below a certain threshold
    threshold = max(anglespec,[],'all')*0.5;  % Set your threshold value here
    local_maxima(anglespec < threshold) = false;
    
    % Step 4: Get indices of remaining maxima
    [row_indices, col_indices] = find(local_maxima);
    
    % Step 5: Map indices to scales and angles
    max_scales = Scales(row_indices);
    max_angles = Angles_pos(col_indices);

    % Step 9: Prepare the list of scale/angle pairs
    peak_list = [max_scales; max_angles];  % Create a list of scale/angle pairs for the upper half

    %% Part 1: Plot coherence with red contour based on the 60% threshold
   
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

        % Define an Angles array corresponding to the indices in the phase_difference 4th dimension
        num_angles = size(phase_difference, 4);  % Number of angles in phase_difference
        Angles = linspace(0, 2*pi, num_angles);  % Create the angle array based on your data structure
        
        % Find the index of the closest angle in the Angles array
        [~, angle_idx] = min(abs(Angles - real_angle));
        
        % Extract the corresponding phase_difference slice for the current scale and angle
        phase_slice = phase_difference(:, :, scale_idx, angle_idx);
     
        % Calculate the mean phase value for the masked region
        phase_mean = mean(phase_slice(coherence_mask), 'omitnan');
        
        % Add the phase mean value to the third row of peak_list
        peak_list(3, i) = phase_mean;
   
    end
    
    %% Calculate Speed from Phase Shift
    
    % Initialize arrays to store results
    num_peaks = size(peak_list, 2);
    speeds = zeros(1, num_peaks);
    distance_shifts = zeros(1, num_peaks);
    
    for i = 1:num_peaks
        % Extract scale and mean phase difference from peak_list
        scale = peak_list(1, i);
        mean_phase_difference = peak_list(3, i);
        
        % Calculate the wavelength in km
        % Since scale is half the wavelength in pixels, and each pixel is pixel_size_km
        wavelength_km = scale * 2 * pixel_size_km;
        
        % Calculate the distance shift in km
        distance_shift_km = mean_phase_difference * wavelength_km / (2 * pi);
        
        % Calculate the speed in km/s
        speed_km_per_s = distance_shift_km / time_interval;
        
        % Store the results
        distance_shifts(i) = distance_shift_km;
        speeds(i) = speed_km_per_s;
    end

    % Add the speeds to the peak_list
    peak_list(4, :) = speeds * 1000; % Convert to m/s
end

function analyze_peak_data(peak_data, num_squares_x, num_squares_y)
    %% Analyze Peak Data Function
    % This function takes the peak_data array and performs statistical analyses
    % and visualizations to identify patterns and correlations in the data.
    %
    % Parameters:
    %   peak_data      - The array containing the peak data.
    %   num_squares_x  - Number of squares along the x-axis.
    %   num_squares_y  - Number of squares along the y-axis.
    %
    % The function generates histograms, scatter plots, heatmaps, PCA analysis,
    % correlation matrices, and time evolution plots.

    % Extract columns from peak_data
    frame_ids    = peak_data(:, 1);
    square_ids   = peak_data(:, 2);
    scales       = peak_data(:, 3);
    angles       = peak_data(:, 4);
    mean_phases  = peak_data(:, 5);
    speeds       = peak_data(:, 6);

    % Remove any rows with NaN values
    valid_idx = ~any(isnan(peak_data), 2);
    frame_ids    = frame_ids(valid_idx);
    square_ids   = square_ids(valid_idx);
    scales       = scales(valid_idx);
    angles       = angles(valid_idx);
    mean_phases  = mean_phases(valid_idx);
    speeds       = speeds(valid_idx);

    %% Histograms of Scales, Angles, and Speeds
    figure;
    histogram(scales, 'BinWidth', 0.5);
    title('Histogram of Scales');
    xlabel('Scale');
    ylabel('Frequency');

    figure;
    histogram(angles, 'BinWidth', pi/18); % Bin width of 10 degrees
    title('Histogram of Angles');
    xlabel('Angle (radians)');
    ylabel('Frequency');

    figure;
    histogram(speeds, 'BinWidth', 0.1);
    title('Histogram of Speeds');
    xlabel('Speed (m/s)');
    ylabel('Frequency');

    %% Scatter Plot of Scales vs Angles Colored by Speed
    figure;
    scatter(angles, scales, 20, speeds, 'filled');
    colorbar;
    title('Scatter Plot of Scales vs Angles');
    xlabel('Angle (radians)');
    ylabel('Scale');
    colormap('jet');

    %% Mean Values per Square
    unique_squares = unique(square_ids);
    num_squares = length(unique_squares);
    mean_scales = zeros(num_squares, 1);
    mean_angles = zeros(num_squares, 1);
    mean_speeds = zeros(num_squares, 1);

    for i = 1:num_squares
        sq_id = unique_squares(i);
        idx = square_ids == sq_id;
        mean_scales(i) = mean(scales(idx));
        mean_angles(i) = mean(angles(idx));
        mean_speeds(i) = mean(speeds(idx));
    end

    %% Heatmaps of Mean Scales, Angles, and Speeds per Square
    mean_scales_matrix = nan(num_squares_y, num_squares_x);
    mean_angles_matrix = nan(num_squares_y, num_squares_x);
    mean_speeds_matrix = nan(num_squares_y, num_squares_x);

    for i = 1:num_squares
        sq_id = unique_squares(i);
        % Map sq_id to grid position
        row = ceil(sq_id / num_squares_x);
        col = mod(sq_id - 1, num_squares_x) + 1;
        mean_scales_matrix(row, col) = mean_scales(i);
        mean_angles_matrix(row, col) = mean_angles(i);
        mean_speeds_matrix(row, col) = mean_speeds(i);
    end

    % Plot heatmap of mean scales
    figure;
    imagesc(mean_scales_matrix);
    colorbar;
    title('Heatmap of Mean Scales per Square');
    xlabel('X Position');
    ylabel('Y Position');
    set(gca, 'YDir', 'normal'); % Correct the Y-axis direction

    % Plot heatmap of mean angles
    figure;
    imagesc(mean_angles_matrix);
    colorbar;
    title('Heatmap of Mean Angles per Square');
    xlabel('X Position');
    ylabel('Y Position');
    set(gca, 'YDir', 'normal');

    % Plot heatmap of mean speeds
    figure;
    imagesc(mean_speeds_matrix);
    colorbar;
    title('Heatmap of Mean Speeds per Square');
    xlabel('X Position');
    ylabel('Y Position');
    set(gca, 'YDir', 'normal');

    %% Principal Component Analysis (PCA)
    % Data matrix: scales, angles, speeds
    X = [scales, angles, speeds];
    [coeff, score, latent, ~, explained] = pca(X);

    % Plot explained variance
    figure;
    pareto(explained);
    xlabel('Principal Component');
    ylabel('Variance Explained (%)');
    title('PCA Explained Variance');

    % Scatter plot of the first two principal components
    figure;
    scatter(score(:, 1), score(:, 2), 20, frame_ids, 'filled');
    colorbar;
    title('PCA - First Two Principal Components');
    xlabel('PC1');
    ylabel('PC2');
    colormap('jet');

    %% Correlation Matrix
    corr_matrix = corrcoef(X);
    figure;
    imagesc(corr_matrix);
    colorbar;
    title('Correlation Matrix');
    xticks(1:3);
    xticklabels({'Scale', 'Angle', 'Speed'});
    yticks(1:3);
    yticklabels({'Scale', 'Angle', 'Speed'});
    set(gca, 'XTickLabelRotation', 45);

    %% Time Evolution of Mean Values
    unique_frames = unique(frame_ids);
    num_frames = length(unique_frames);
    mean_scales_frame = zeros(num_frames, 1);
    mean_angles_frame = zeros(num_frames, 1);
    mean_speeds_frame = zeros(num_frames, 1);

    for i = 1:num_frames
        fr_id = unique_frames(i);
        idx = frame_ids == fr_id;
        mean_scales_frame(i) = mean(scales(idx));
        mean_angles_frame(i) = mean(angles(idx));
        mean_speeds_frame(i) = mean(speeds(idx));
    end

    % Plot mean scale over time
    figure;
    plot(unique_frames, mean_scales_frame, 'o-');
    title('Mean Scale over Frames');
    xlabel('Frame ID');
    ylabel('Mean Scale');

    % Plot mean angle over time
    figure;
    plot(unique_frames, mean_angles_frame, 'o-');
    title('Mean Angle over Frames');
    xlabel('Frame ID');
    ylabel('Mean Angle (radians)');

    % Plot mean speed over time
    figure;
    plot(unique_frames, mean_speeds_frame, 'o-');
    title('Mean Speed over Frames');
    xlabel('Frame ID');
    ylabel('Mean Speed (m/s)');

    %% Correlation Analysis Between Variables
    % Correlation between scales and speeds
    [r_scale_speed, p_scale_speed] = corr(scales, speeds, 'Type', 'Pearson');
    fprintf('Correlation between scales and speeds: r = %.2f, p = %.4f\n', r_scale_speed, p_scale_speed);

    % Correlation between angles and speeds
    [r_angle_speed, p_angle_speed] = corr(angles, speeds, 'Type', 'Pearson');
    fprintf('Correlation between angles and speeds: r = %.2f, p = %.4f\n', r_angle_speed, p_angle_speed);

    % Correlation between scales and angles
    [r_scale_angle, p_scale_angle] = corr(scales, angles, 'Type', 'Pearson');
    fprintf('Correlation between scales and angles: r = %.2f, p = %.4f\n', r_scale_angle, p_scale_angle);

    %% 3D Scatter Plot of Scales, Angles, and Speeds
    figure;
    scatter3(scales, angles, speeds, 20, frame_ids, 'filled');
    colorbar;
    title('3D Scatter Plot of Scales, Angles, and Speeds');
    xlabel('Scale');
    ylabel('Angle (radians)');
    zlabel('Speed (m/s)');
    colormap('jet');
    grid on;

    %% Additional Analysis
    % You can extend this function with more specific analyses, such as:
    % - Investigating specific squares with unusual patterns.
    % - Analyzing the distribution of mean phases.
    % - Performing cluster analysis to identify groups of similar data points.
    % - Examining the time evolution of patterns in specific squares.

    disp('Analysis complete.');
end

function plot_waverose(frame_id, square_id)
    %% Plot Waverose Function
    % This function generates an advanced rose plot for a specific frame and square.
    % It overlays two half circles: one representing power and the other representing coherence.
    
    % Parameters:
    %   frame_id   - The ID of the frame to analyze.
    %   square_id  - The ID of the square to analyze.
    
    % Ensure that the necessary variables and functions are accessible
    % If they are not in the workspace, you may need to pass them as additional arguments
    % or define them within this function.
    
    % Load necessary global parameters or define them here
    % For example, video_file, pixel_size_km, Scales, Angles, etc.
    % Assuming these variables are accessible in the workspace
    
    %% Load Parameters from Base Workspace
    % Retrieve variables from the base workspace
    video_file = evalin('base', 'video_file');
    shrinkfactor = evalin('base', 'shrinkfactor');
    invshrinkfactor = evalin('base', 'invshrinkfactor');
    pixel_size_km = evalin('base', 'pixel_size_km');
    window_function = evalin('base', 'window_function');
    Scales = evalin('base', 'Scales');
    Angles = evalin('base', 'Angles');
    time_interval = evalin('base', 'time_interval');
    squares = evalin('base', 'squares');
    Preprocess_Flag = evalin('base', 'Preprocess_Flag');
    
    %% Read Video Frames
    % Initialize video reader
    v = VideoReader(video_file);
    
    % Check if frame_id is valid
    if frame_id < 1 || frame_id >= v.NumFrames
        error('Invalid frame_id. Must be between 1 and %d.', v.NumFrames - 1);
    end
    
    % Read the specified frames
    frame1 = read(v, frame_id);
    frame2 = read(v, frame_id + 1);
    
    % Resize frames according to shrink factor
    frame1 = imresize(frame1, invshrinkfactor);
    frame2 = imresize(frame2, invshrinkfactor);
    
    % Convert to grayscale if necessary and extract red channel
    if size(frame1, 3) > 1
        frame1 = double(frame1(:, :, 1));
        frame2 = double(frame2(:, :, 1));
    else
        frame1 = double(frame1);
        frame2 = double(frame2);
    end
    
    % Preprocess frames if needed
    if Preprocess_Flag
        frame1 = preprocess_img(frame1);
        frame2 = preprocess_img(frame2);
    end
    
    %% Extract the Specified Square
    % Find the square with the given square_id
    idx = find([squares.index] == square_id);
    if isempty(idx)
        error('Invalid square_id. Square not found.');
    end
    square = squares(idx);
    
    % Extract the square from both frames
    img1 = frame1(square.y_range, square.x_range);
    img2 = frame2(square.y_range, square.x_range);
    
    % Apply windowing to reduce edge effects
    img1_windowed = apply_window(img1, window_function);
    img2_windowed = apply_window(img2, window_function);
    
    %% Perform Wavelet Transform
    % Perform 2D Continuous Wavelet Transform on both frames
    cwt1 = cwtft2(img1_windowed, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
    spec1 = squeeze(cwt1.cfs);
    
    cwt2 = cwtft2(img2_windowed, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
    spec2 = squeeze(cwt2.cfs);
    
    % Compute the cross-wavelet spectrum (XWT)
    xwt = spec1 .* conj(spec2);
    
    % Compute coherence and phase difference
    power1 = abs(spec1).^2;
    power2 = abs(spec2).^2;
    coherence = abs(xwt);
    phase_difference = angle(xwt);
    
    %% Find Peaks and Speeds
    % Use the existing function to find peaks and calculate speeds
    peak_list = find_peaks_and_speeds(coherence, phase_difference, Scales, Angles, pixel_size_km, time_interval);
    
    %% Plot Advanced Rose Plot
    % Generate the advanced rose plot with power and coherence
    % Overlay the peaks on the plot
    
    % Calculate inner power and coherence for plotting
    buffer = round(max(Scales));
    innerpower = squeeze(mean(mean(power1(buffer:end-buffer, buffer:end-buffer, :, :))));
    innercoherence = squeeze(mean(mean(coherence(buffer:end-buffer, buffer:end-buffer, :, :))));
    
    % Normalize power and coherence by scale
    mean_power_byscale = mean(innerpower, 2);
    mean_coherence_byscale = mean(innercoherence, 2);
    
    anglespec_power = innerpower ./ mean_power_byscale;
    anglespec_coherence = innercoherence ./ mean_coherence_byscale;
    
    % Define angles for upper and lower halves
    Angles_pos = Angles;
    Angles_neg = Angles + pi;
    
    % Prepare data for plotting
    [Theta_pos, R_pos] = meshgrid(Angles_pos, Scales);
    [X_pos, Y_pos] = pol2cart(Theta_pos, R_pos);
    
    [Theta_neg, R_neg] = meshgrid(Angles_neg, Scales);
    [X_neg, Y_neg] = pol2cart(Theta_neg, R_neg);
    
    % Plot the advanced rose plot
    figure;
    % Plot power (upper half)
    ax1 = axes;
    pcolor(ax1, X_pos, Y_pos, anglespec_power);
    shading interp;
    colormap(ax1, 'parula');
    axis equal;
    set(ax1, 'Position', [0.1, 0.1, 0.75, 0.75]);
    ax1.XTick = [];
    ax1.YTick = [];
    hold on;
    
    % Plot coherence (lower half)
    ax2 = axes;
    pcolor(ax2, X_neg, Y_neg, anglespec_coherence);
    shading interp;
    colormap(ax2, 'hot');
    axis equal;
    set(ax2, 'Position', [0.1, 0.1, 0.75, 0.75]);
    ax2.XTick = [];
    ax2.YTick = [];
    set(ax2, 'Color', 'none');
    linkaxes([ax1, ax2]);
    hold on;
    
    % Overlay peaks on the coherence plot
    if ~isempty(peak_list)
        max_scales = peak_list(1, :);
        max_angles = peak_list(2, :) + pi;  % Adjust angles for lower half
        [peak_X, peak_Y] = pol2cart(max_angles, max_scales);
        plot(ax2, peak_X, peak_Y, 'k*', 'MarkerSize', 10);
        % Annotate peaks with speeds
        for i = 1:length(peak_X)
            text(ax2, peak_X(i) * 1.05, peak_Y(i) * 1.05, sprintf('%.2f m/s', peak_list(4, i)), 'Color', 'k', 'FontSize', 10);
        end
    end
    
    % Adjust axes limits
    xlim(ax1, [min(X_pos(:)) - 1, max(X_pos(:)) + 1]);
    ylim(ax1, [min(Y_neg(:)) - 1, max(Y_pos(:)) + 1]);
    
    %% Add Radial Rings and Angle Labels
    % Add radial rings corresponding to scales
    ring_radii = Scales;
    for i = 1:length(ring_radii)
        theta_ring = linspace(0, 2 * pi, 100);
        [x_ring, y_ring] = pol2cart(theta_ring, ring_radii(i));
        plot(ax1, x_ring, y_ring, 'k--');
        plot(ax2, x_ring, y_ring, 'k--');
        % Add scale labels
        text(ax1, ring_radii(i) * 1.05, 0, num2str(ring_radii(i), '%.2f'), 'HorizontalAlignment', 'left');
    end
    
    % Add angle lines and labels
    angle_ticks = linspace(0, 2 * pi, 13);  % Every 30 degrees
    angle_labels = {'0', '\pi/6', '\pi/3', '\pi/2', '2\pi/3', '5\pi/6', '\pi', '7\pi/6', '4\pi/3', '3\pi/2', '5\pi/3', '11\pi/6', '2\pi'};
    for i = 1:length(angle_ticks)
        angle_rad = angle_ticks(i);
        x_line = [0, max(Scales) * cos(angle_rad)];
        y_line = [0, max(Scales) * sin(angle_rad)];
        plot(ax1, x_line, y_line, 'k--');
        plot(ax2, x_line, y_line, 'k--');
        text(ax1, x_line(2) * 1.1, y_line(2) * 1.1, angle_labels{i}, 'HorizontalAlignment', 'center');
    end
    
    %% Add Colorbars
    original_pos = get(ax1, 'Position');
    c1 = colorbar(ax1, 'eastoutside');
    c1_pos = get(c1, 'Position');
    c1_pos(1) = c1_pos(1) + 0.05;
    set(c1, 'Position', c1_pos);
    set(ax1, 'Position', original_pos);
    ylabel(c1, 'Power');
    
    original_pos = get(ax2, 'Position');
    c2 = colorbar(ax2, 'westoutside');
    c2_pos = get(c2, 'Position');
    c2_pos(1) = c2_pos(1) - 0.05;
    set(c2, 'Position', c2_pos);
    set(ax2, 'Position', original_pos);
    ylabel(c2, 'Coherence');
    
    %% Set Titles
    sgtitle(sprintf('Advanced Rose Plot for Frame %d and Square %d', frame_id, square_id));
    
end