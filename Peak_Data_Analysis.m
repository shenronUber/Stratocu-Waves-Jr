% --- Code to analyze data and visualize displacements ---

% 1. Analysis of average coherence
coherence_values = peak_data(:, 7);

% Plot histogram of coherence values with a logarithmic y-axis
figure;
histogram(coherence_values, 'Normalization', 'probability'); % Normalize for better log visualization
set(gca, 'YScale', 'log'); % Set the y-axis to logarithmic scale
xlabel('Coherence Value');
ylabel('Frequency (log)');
title('Histogram of Coherence Values (Logarithmic Scale)');

% 2. Correlation between wavelength and coherence
wavelengths = peak_data(:, 8);
coherence_values = peak_data(:, 7); % Ensure this contains coherence values

figure;
boxplot(coherence_values, wavelengths, 'PlotStyle', 'compact', 'Symbol', '.');
xlabel('Wavelength (km)');
ylabel('Coherence Value');
title('Coherence Distribution by Wavelength');

% 3. Histogram of coherence values normalized by the mean for each wavelength

% Calculate coherence normalized by mean for each wavelength
unique_wavelengths = unique(wavelengths);
normalized_coherence_values = zeros(size(coherence_values));

for i = 1:length(unique_wavelengths)
    current_wavelength = unique_wavelengths(i);
    indices = wavelengths == current_wavelength;
    mean_coherence = mean(coherence_values(indices));
    normalized_coherence_values(indices) = coherence_values(indices) / mean_coherence;
end

% 4. Correlation between wavelength and normalized coherence
figure;
boxplot(normalized_coherence_values, wavelengths, 'PlotStyle', 'compact', 'Symbol', '.');
xlabel('Wavelength (km)');
ylabel('Normalized Coherence');
title('Normalized Coherence Distribution by Wavelength');

% 5. Filter coherence values based on outliers for each wavelength

% Calculate quartiles and outliers for each unique wavelength
filtered_coherence_values = normalized_coherence_values; % Copy of original values

for i = 1:length(unique_wavelengths)
    current_wavelength = unique_wavelengths(i);
    indices = wavelengths == current_wavelength;
    
    % Compute quartiles and IQR for this wavelength
    Q1 = prctile(normalized_coherence_values(indices), 25);
    Q3 = prctile(normalized_coherence_values(indices), 75);
    IQR = Q3 - Q1;
    
    % Detect and exclude outliers
    outliers = normalized_coherence_values > (Q3 + 1.5 * IQR) | normalized_coherence_values < (Q1 - 1.5 * IQR);
    filtered_coherence_values(indices & outliers) = NaN; % Replace outliers with NaN
end


%% Parameters and Setup

% Start timing the entire process
total_time = tic;

% Directory containing NetCDF files
data_dir = 'downloaded_data\VIS_2023_11_01-07';

% Number of frames to process
% Set to Inf to process all frames, or specify a number, e.g., 10 for the first 10 frames
num_frames = 2;

% Get a list of NetCDF files in the directory
nc_files = dir(fullfile(data_dir, '*.nc'));

% Adjust num_frames based on available files
if num_frames == Inf
    num_frames = numel(nc_files); % Process all files
else
    num_frames = min(num_frames, numel(nc_files)); % Process up to the specified limit
end

% Ensure there are enough files for processing
if num_frames < 2
    error('At least two NetCDF files are required for processing.');
end

% Shrinking factor for resizing
shrinkfactor = 2;
invshrinkfactor = 1 / shrinkfactor;

% Pixel size in kilometers
degrees_per_pixel = 0.04;
km_per_degree = 111.32; % Average km per degree on Earth's surface
original_pixel_size_km = degrees_per_pixel * km_per_degree; % Resulting pixel size in km (4.4528 km)
pixel_size_km = original_pixel_size_km * shrinkfactor; % Adjusted pixel size after shrinking

% Square size (5x5 degrees)
square_size_deg = 5; 
square_size_km = square_size_deg * km_per_degree; % Total km per square
square_size_px = round(square_size_km / pixel_size_km);

% Brightness thresholds for filtering
brightness_threshold = 0.00001; % Minimum mean brightness for inclusion
std_threshold = 10;      % Maximum standard deviation for inclusion

% Wavelet transform parameters
NANGLES = 24;
Angles = 0:pi/(NANGLES-1):pi;

% Define scales in km (10 km to 500 km range, divided by 2 for half-wavelength scales)
min_scale_km = 10/2;
max_scale_km = 500/2;
NSCALES = 20;
Scales_km = logspace(log10(min_scale_km), log10(max_scale_km), NSCALES);

% Convert scales to pixels
Scales = Scales_km / pixel_size_km;

% Windowing parameters for edge effects
window_buffer = 10; % Buffer size in pixels around frame edges

% Preprocessing toggle
Preprocess_Flag = 1; % 1 for on, 0 for off

% Time interval between frames (adjust as necessary)
time_interval = 1800; % Assuming 30 minutes in seconds

% Metadata saving toggle
Save_Metadata_Flag = 0; % 1 to save metadata, 0 to skip

%% Read Data and Initialize

% Read the first file to get dimensions
first_file_path = fullfile(data_dir, nc_files(1).name);
frame1 = readVariableFromFile(first_file_path);

% Transpose the frame for portrait orientation
frame1 = frame1';

% Resize the frame using the shrinking factor
if shrinkfactor ~= 1
    frame1 = imresize(frame1, invshrinkfactor);
end

[frame_height, frame_width] = size(frame1);

% Adjust frame dimensions to exclude the buffer
x_buffer_range = (window_buffer + 1):(frame_width - window_buffer);
y_buffer_range = (window_buffer + 1):(frame_height - window_buffer);

adjusted_frame_width = length(x_buffer_range);
adjusted_frame_height = length(y_buffer_range);

% Calculate number of squares along width and height
num_squares_x = ceil(adjusted_frame_width / square_size_px);
num_squares_y = ceil(adjusted_frame_height / square_size_px);

% Generate Squares Coordinates

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
if Preprocess_Flag == 1
    frame1 = preprocess_img(frame1);
end

% Define the radial window parameters
radius_factor = 1.1;  % 110% of the image radius
decay_rate = 0.05;    % Decay rate for windowing

% Apply the radial window
frame1_windowed = apply_radial_window(frame1, radius_factor, decay_rate);

imshow(frame1_windowed, [0, 1]);
title('Windowed Frame 1');
colormap('gray');

hold on;
for idx = 1:length(squares)
    rectangle('Position', [squares(idx).x_range(1), squares(idx).y_range(1), ...
        length(squares(idx).x_range), length(squares(idx).y_range)], ...
        'EdgeColor', 'r', 'LineWidth', 1);
    text(squares(idx).x_range(1), squares(idx).y_range(1) - 5, ...
        sprintf('%d', squares(idx).index), 'Color', 'yellow', 'FontSize', 8);
end
title('Frame 1 with Squares');
hold off;

%%
% 5. Filtering the data

% Assume filtered_peak_data is a copy of peak_data
filtered_peak_data = peak_data;

% Define the desired wavelength in km
desired_wavelength = 50; % For example, 100 km

% Define the tolerance in km (e.g., ±10% of the desired wavelength)
wavelength_tolerance = desired_wavelength * 0.1;

% Index of the wavelength column in peak_data
wavelength_column = 8; % Ensure the wavelength is in column 8

% Compute the absolute difference between peak wavelengths and the desired wavelength
wavelength_diff = abs(filtered_peak_data(:, wavelength_column) - desired_wavelength);

% Find indices of peaks where the difference is less than or equal to the tolerance
valid_indices = wavelength_diff <= wavelength_tolerance;

% Create a new table with peaks filtered by wavelength
desired_wavelength_peaks = filtered_peak_data(valid_indices, :);

% Threshold for coherence (modify as needed)
coherence_threshold = 0 * 0.5 * mean(filtered_peak_data(:, 7));

% Replace coherence values in filtered_peak_data with filtered_coherence_values
filtered_peak_data(:, 7) = filtered_coherence_values;

% Find valid indices where coherence values are not NaN
valid_indices = ~isnan(filtered_coherence_values);

% Filter the data to keep only valid rows
filtered_peak_data = filtered_peak_data(valid_indices, :);

% Plot the histogram of filtered coherence values to check outlier reduction
figure;
histogram(filtered_coherence_values(valid_indices), 'Normalization', 'probability');
set(gca, 'YScale', 'log');
xlabel('Coherence Value (after outlier filtering)');
ylabel('Frequency (log)');
title('Histogram of Filtered Coherence Values (Logarithmic Scale)');

% 4. Tracking peaks between frames
num_frames = max(peak_data(:, 1));
num_squares = max(peak_data(:, 2));

% Initialize peaks_by_square
peaks_by_square = cell(num_squares, 1);
for s = 1:num_squares
    peaks_by_square{s} = cell(num_frames, 1);
end

% Populate peaks_by_square with peaks filtered by wavelength
for i = 1:size(desired_wavelength_peaks, 1)
    frame_idx = desired_wavelength_peaks(i, 1);
    square_idx = desired_wavelength_peaks(i, 2);
    peaks_by_square{square_idx}{frame_idx} = desired_wavelength_peaks(i, :);
end

% Calculate displacements between frames
displacements = [];

for s = 1:num_squares
    square_peaks = peaks_by_square{s};
    square_displacements = [];
    
    for frame_idx = 1:(num_frames - 1)
        peak_current = square_peaks{frame_idx};
        peak_next = square_peaks{frame_idx + 1};
        
        if ~isempty(peak_current) && ~isempty(peak_next)
            angle_current = peak_current(4);
            angle_next = peak_next(4);
            angle_diff = angle_next - angle_current;
            angle_diff = atan2(sin(angle_diff), cos(angle_diff));
            
            wavelength_current = peak_current(8);
            wavelength_next = peak_next(8);
            wavelength_diff = abs(wavelength_next - wavelength_current);
            
            displacement = struct();
            displacement.square_idx = s;
            displacement.frame_idx = frame_idx;
            displacement.angle_current = angle_current;
            displacement.angle_next = angle_next;
            displacement.angle_diff = angle_diff;
            displacement.wavelength_current = wavelength_current;
            displacement.wavelength_next = wavelength_next;
            displacement.wavelength_diff = wavelength_diff;
            displacement.coherence_current = peak_current(7);
            displacement.coherence_next = peak_next(7);
            
            square_displacements = [square_displacements; displacement];
        end
    end
    displacements = [displacements; square_displacements];
end

% 5. Filter significant displacements
angle_threshold = pi / 8;
wavelength_threshold = 0.1 * mean(wavelengths);

significant_displacements = [];
for i = 1:length(displacements)
    displacement = displacements(i);
    if abs(displacement.angle_diff) <= angle_threshold && displacement.wavelength_diff <= wavelength_threshold
        significant_displacements = [significant_displacements; displacement];
    end
end

% 6. Visualize arrows for the desired wavelength

% Load and preprocess frame1_display as before
frame1_path = fullfile(data_dir, nc_files(1).name);
frame1 = readVariableFromFile(frame1_path);
frame1 = frame1';

if shrinkfactor ~= 1
    frame1 = imresize(frame1, invshrinkfactor);
end

if Preprocess_Flag
    frame1_display = preprocess_img(frame1);
else
    frame1_display = frame1;
end

frame1_display = apply_radial_window(frame1_display, radius_factor, decay_rate);

figure;
imshow(frame1_display, []);
hold on;

% Threshold for visualizing coherence
coherence_display_threshold = coherence_threshold;

for s = 1:num_squares
    peak = peaks_by_square{s}{1}; % Now contains peaks for the desired wavelength
    if ~isempty(peak) && peak(7) >= coherence_display_threshold
        % Extract peak parameters
        angle = peak(4);
        speed = peak(6);
        wavelength = peak(8);

        % Calculate arrow properties
        square = squares(s);
        x_center = mean(square.x_range);
        y_center = mean(square.y_range);

        arrow_length = 20;
        dx = arrow_length * cos(angle);
        dy = arrow_length * sin(angle);

        quiver(x_center, y_center, dx, dy, 'LineWidth', 2, 'Color', 'red', 'MaxHeadSize', 2);

        text(x_center + dx, y_center + dy, sprintf('V: %.1f m/s\nλ: %.1f km', speed, wavelength), ...
            'Color', 'yellow', 'FontSize', 8);
    end
end

title(sprintf('Displacement Directions Detected for λ ≈ %.1f km', desired_wavelength));
hold off;


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
