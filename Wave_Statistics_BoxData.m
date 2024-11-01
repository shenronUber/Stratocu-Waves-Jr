%% Parameters and Setup

% Start timing the entire process
total_time = tic;

% Directory containing NetCDF files
data_dir = 'downloaded_data\VIS_2023_11_01-07';

% Specify the number of frames to process
% Set to Inf to process all frames, or specify a number, e.g., 10 for the first 10 frames
num_frames = 2;

% Get a list of NetCDF files in the directory
nc_files = dir(fullfile(data_dir, '*.nc'));

% Adjust num_frames based on available files
if num_frames == Inf
    num_frames = numel(nc_files); % Process all frames
else
    num_frames = min(num_frames, numel(nc_files)); % Process up to the specified limit
end

% Check if there are enough files to process
if num_frames < 2
    error('At least two NetCDF files are required for processing.');
end

% Shrinking factor
shrinkfactor = 2;
invshrinkfactor = 1 / shrinkfactor;

% Dynamic pixel size (km per pixel)
degrees_per_pixel = 0.04;
km_per_degree = 111.32; % Average km per degree on Earth's surface
original_pixel_size_km = degrees_per_pixel * km_per_degree; % Resulting pixel size in km (4.4528 km)
pixel_size_km = original_pixel_size_km * shrinkfactor; % Adjusted pixel size due to shrinking

% Square size in degrees (does not depend on shrinkfactor)
square_size_deg = 5; % 5x5 degrees squares

% Calculate square size in km (independent of shrinkfactor)
square_size_km = square_size_deg * km_per_degree; % Total km per square

% Calculate square size in pixels based on original pixel size
square_size_px = round(square_size_km / pixel_size_km);

% Brightness thresholds (adjust as needed)
brightness_threshold = 0.00001; % Mean brightness below which squares are ignored
std_threshold = 10;      % Standard deviation above which squares are ignored

% Wavelet transform parameters
NANGLES = 24;
Angles = 0:pi/(NANGLES-1):pi;

% Define scales in km (10 km to 100 km range)
min_scale_km = 10;
max_scale_km = 500;
NSCALES=20;
Scales_km = logspace(log10(min_scale_km), log10(max_scale_km), NSCALES);

% Convert scales to pixels (independent of shrinkfactor)
Scales = Scales_km / pixel_size_km;

% Windowing parameters (independent of shrinkfactor)
window_buffer = 10; % Buffer size in pixels around the frame edges

% Preprocessing flag
Preprocess_Flag = 1; % 1 is on / 0 is off

% Time interval between frames (adjust if necessary)
time_interval = 1800; % Assuming 30 minutes in seconds

% Metadata saving flag
Save_Metadata_Flag = 1; % 1 to save metadata, 0 to skip

%% Read Data and Initialize

% Read the first file to get dimensions
first_file_path = fullfile(data_dir, nc_files(1).name);
frame1 = readVariableFromFile(first_file_path);

% Transpose the frame to switch x and y dimensions (portrait orientation)
frame1 = frame1';

% Resize the frame according to the shrink factor if needed
if shrinkfactor ~= 1
    frame1 = imresize(frame1, invshrinkfactor);
end

[frame_height, frame_width] = size(frame1);

% Adjust frame dimensions by removing the buffer
x_buffer_range = (window_buffer + 1):(frame_width - window_buffer);
y_buffer_range = (window_buffer + 1):(frame_height - window_buffer);

adjusted_frame_width = length(x_buffer_range);
adjusted_frame_height = length(y_buffer_range);

% Calculate number of squares along width and height
num_squares_x = ceil(adjusted_frame_width / square_size_px);
num_squares_y = ceil(adjusted_frame_height / square_size_px);

% Generate Squares Coordinates

% Generate coordinates for squares
squares = [];
idx = 1;
for i = 1:num_squares_y
    for j = 1:num_squares_x
        % Calculate pixel indices
        x_start = floor((j - 1) * adjusted_frame_width / num_squares_x) + 1;
        y_start = floor((i - 1) * adjusted_frame_height / num_squares_y) + 1;
        x_end = floor(j * adjusted_frame_width / num_squares_x);
        y_end = floor(i * adjusted_frame_height / num_squares_y);

        % Ensure valid ranges
        if x_end > x_start && y_end > y_start
            squares(idx).x_range = x_buffer_range(x_start:x_end);
            squares(idx).y_range = y_buffer_range(y_start:y_end);
            squares(idx).index = idx;
            idx = idx + 1;
        end
    end
end

% Display Frame 1 with Squares

figure;
% Preprocess frames if needed
if Preprocess_Flag==1
    frame1 = preprocess_img(frame1);
end

% Define the parameters of the radial window
radius_factor = 1.1;  % 110% of the image radius (to keep the whole image)
decay_rate = 0.05;    % Control parameter for decay rate

% Apply windowing only at the edges
frame1_windowed = apply_radial_window(frame1, radius_factor, decay_rate);

imshow(frame1_windowed, [0,1]);
title('Windowed Frame 1');
colormap('gray');

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

% Adjusted number of entries based on the number of scales
total_entries = (num_frames - 1) * length(squares) * NSCALES;

% Initialize peak_data: columns are [frame_idx, square_idx, scale, angle, phase_difference, speed, coherence_value, wavelength_km]
peak_data = NaN(total_entries, 8); % Now 8 columns to include wavelength

% Initialize a counter for the number of entries
entry_counter = 0;

% Initialize a structure array to store metadata if needed
if Save_Metadata_Flag
    square_metadata = struct('frame_idx', {}, 'square_idx', {}, 'spec1', {});
end

% Loop over files (frames)
for frame_idx = 1:num_frames - 1
    fprintf('Processing frame %d/%d...\n', frame_idx, num_frames);

    % Measure time for processing each frame (optional)
    frame_time = tic; % Start timing for a specific frame

    % Read data from NetCDF files
    file1_path = fullfile(data_dir, nc_files(frame_idx).name);
    file2_path = fullfile(data_dir, nc_files(frame_idx + 1).name);

    frame1 = readVariableFromFile(file1_path);
    frame2 = readVariableFromFile(file2_path);

    % Transpose the frames to switch x and y dimensions (portrait orientation)
    frame1 = frame1';
    frame2 = frame2';

    % Resize frames according to shrink factor if needed
    if shrinkfactor ~= 1
        frame1 = imresize(frame1, invshrinkfactor);
        frame2 = imresize(frame2, invshrinkfactor);
    end

    % Preprocess frames if needed
    if Preprocess_Flag
        frame1 = preprocess_img(frame1);
        frame2 = preprocess_img(frame2);
    end

    % Apply windowing to the entire frame
    frame1_windowed = apply_radial_window(frame1, radius_factor, decay_rate);
    frame2_windowed = apply_radial_window(frame2, radius_factor, decay_rate);

    % Adjust frame dimensions after applying buffer
    frame1_windowed = frame1_windowed(y_buffer_range, x_buffer_range);
    frame2_windowed = frame2_windowed(y_buffer_range, x_buffer_range);

    % Perform wavelet transform on the entire frames
    cwt1_full = cwtft2(frame1_windowed, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
    spec1_full = squeeze(cwt1_full.cfs);

    cwt2_full = cwtft2(frame2_windowed, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
    spec2_full = squeeze(cwt2_full.cfs);

    % Loop over squares
    for s_idx = 1:length(squares)
        square = squares(s_idx);

        % Extract the square from both frames' wavelet coefficients
        % Get the indices for the square in the cwt coefficient arrays
        x_range = square.x_range - window_buffer;
        y_range = square.y_range - window_buffer;

        % Extract the cwt coefficients for the square
        spec1 = spec1_full(y_range, x_range, :, :);
        spec2 = spec2_full(y_range, x_range, :, :);

        % Optionally save the square's metadata (CWT coefficients from spec1)
        if Save_Metadata_Flag == 1
            % Save the CWT coefficients for this square (only from spec1)
            metadata_entry.frame_idx = frame_idx;
            metadata_entry.square_idx = s_idx;
            metadata_entry.spec1 = spec1; % Only saving spec1 (coefficients of the first frame)
            square_metadata = [square_metadata; metadata_entry];
        end

        % Check brightness and standard deviation on the original frames
        img1 = frame1_windowed(y_range, x_range);
        mean_brightness = mean(img1(:));
        std_brightness = std(img1(:));

        if mean_brightness < brightness_threshold || std_brightness > std_threshold
            % Ignore this square
            continue;
        end

        % Proceed with the analysis using the extracted CWT coefficients

        % Compute the cross-wavelet spectrum (XWT)
        xwt = spec1 .* conj(spec2);

        % Compute coherence and phase difference
        power1 = abs(spec1).^2;
        power2 = abs(spec2).^2;
        coherence = abs(xwt);
        phase_difference = angle(xwt);

        % Find peaks in power spectrum and calculate speeds
        peak_list = find_peaks_and_speeds(coherence, phase_difference, Scales, Angles, pixel_size_km, time_interval);

        % Number of scales
        n_scales = size(peak_list, 1);

        % Prepare data to store
        data_to_store = [repmat(frame_idx, n_scales, 1), repmat(s_idx, n_scales, 1), peak_list];

        % Update entry counter
        idx_start = entry_counter + 1;
        idx_end = entry_counter + n_scales;

        % Ensure we do not exceed preallocated array size
        if idx_end > total_entries
            % Expand the array if necessary
            extra_entries = total_entries; % Double the size
            peak_data = [peak_data; NaN(extra_entries, 8)]; % Now 8 columns to match new structure
            total_entries = total_entries + extra_entries;
        end

        % Store data
        peak_data(idx_start:idx_end, :) = data_to_store;

        % Update the entry counter
        entry_counter = idx_end;
    end

    % Display time taken to process this frame
    fprintf('Time to process frame %d: %.2f seconds.\n', frame_idx, toc(frame_time));
end

% Remove unused (NaN) entries from peak_data
peak_data = peak_data(1:entry_counter, :);

% Display total execution time
fprintf('Total execution time: %.2f seconds.\n', toc(total_time));

%% Analyze the data

% to be done

%% Supporting Functions

function data = readVariableFromFile(filePath)
    % Read the appropriate variable ('CMI' or 'Rad') from a NetCDF file
    info = ncinfo(filePath);
    variables = {info.Variables.Name};

    if any(strcmp(variables, 'CMI'))
        variableName = 'CMI';
    elseif any(strcmp(variables, 'Rad'))
        variableName = 'Rad';
    else
        % List all available variables in the file for reference
        warning('No suitable variable found in file %s. Available variables: %s.', ...
                filePath, strjoin(variables, ', '));
        error('No suitable variable to read.');
    end

    % Read the selected data
    data = ncread(filePath, variableName);
end

function img_processed = preprocess_img(img)
    % Preprocess the image by truncating intensity values to reduce contrast.
    
    % Calculate the 10th and 99th percentiles of the image intensity.
    lower_bound = prctile(img(:), 10);
    upper_bound = prctile(img(:), 99);
    
    % Truncate values below the 20th percentile to the lower bound
    % and values above the 80th percentile to the upper bound.
    img_processed = img;
    img_processed(img < lower_bound) = lower_bound;
    img_processed(img > upper_bound) = upper_bound;
    
    % Normalize to the [0, 1] range.
    img_processed = (img_processed - lower_bound) / (upper_bound - lower_bound);
end

function peak_list = find_peaks_and_speeds(coherence, phase_difference, Scales, Angles, pixel_size_km, time_interval)
    % Initialize peak_list to collect data for each scale
    num_scales = length(Scales);
    peak_list = zeros(num_scales, 6); % Columns: [Scale; Angle; MeanPhaseDiff; Speed_m_per_s; CoherenceValue; Wavelength_km]

    % Loop over each scale
    for scale_idx = 1:num_scales
        % Extract coherence at the current scale across all angles
        % Average coherence over x and y dimensions for each angle
        coherence_scale = squeeze(mean(mean(coherence(:, :, scale_idx, :), 1, 'omitnan'), 2, 'omitnan'));  % (num_angles x 1)

        % Find the angle index with the maximum coherence value at this scale
        [max_coherence_value, angle_idx] = max(coherence_scale);

        % Get the corresponding angle
        angle = Angles(angle_idx);

        % Extract the phase difference and coherence slices at this scale and angle
        phase_slice = phase_difference(:, :, scale_idx, angle_idx);
        coherence_slice = coherence(:, :, scale_idx, angle_idx);

        % Compute the mean phase difference over the square
        % Optionally, apply a coherence threshold to mask out low coherence regions
        coherence_threshold = 0.5 * max(coherence_slice(:));
        coherence_mask = coherence_slice >= coherence_threshold;
        if any(coherence_mask(:))
            mean_phase_diff = mean(phase_slice(coherence_mask), 'omitnan');
        else
            mean_phase_diff = NaN; % No significant coherence regions
        end

        % Calculate speed
        wavelength_km = Scales(scale_idx) * 2 * pixel_size_km;  % Wavelength in km
        distance_shift_km = (mean_phase_diff * wavelength_km) / (2 * pi);
        speed_m_per_s = (distance_shift_km / time_interval) * 1000;  % Convert km/s to m/s

        % Store data
        peak_list(scale_idx, :) = [Scales(scale_idx), angle, mean_phase_diff, speed_m_per_s, max_coherence_value, wavelength_km];
    end
end

function img_windowed = apply_radial_window(img, radius_factor, decay_rate)
    % Apply a radial windowing effect that attenuates the image from the edges towards the center
    % Parameters:
    % - radius_factor: controls how quickly the window decays from the center (e.g., 1.1 means 110% of the image radius)
    % - decay_rate: controls the steepness of the attenuation (larger value makes the transition steeper)

    [rows, cols] = size(img);

    % Compute the center of the image
    center_x = cols / 2;
    center_y = rows / 2;

    % Create a meshgrid to calculate the distance from the center
    [x, y] = meshgrid(1:cols, 1:rows);

    % Calculate the radial distance from the center for each pixel
    distances = sqrt((x - center_x).^2 + (y - center_y).^2);

    % Calculate the maximum distance from the center (i.e., the radius of the window)
    max_distance = radius_factor * min(center_x, center_y);  % Factor of the image size

    % Create the radial window using a smooth logistic function
    window = 1 ./ (1 + exp(decay_rate * (distances - max_distance)));

    % Apply the window to the image
    img_windowed = img .* window;
end
