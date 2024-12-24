%% Parameters and Setup (First Part)

startDate = datetime(2023, 10, 12, 1, 0, 0);
endDate = datetime(2023, 10, 12, 9, 0,0 );
downloadDir = 'downloaded_data/IR';
methodName = 'highpass_50'; % e.g., 'raw_normalized' or 'highpass_20_sqrt'

%% List and Select Files Within Date Range
ncFiles = dir(fullfile(downloadDir, '*.nc*'));  % Includes .nc and .nc.4km files
if isempty(ncFiles)
    error('No .nc files found in the download directory.');
end

fileNames = {};
fileTimestamps = datetime([], 'ConvertFrom', 'datenum');

for i = 1:length(ncFiles)
    fileName = ncFiles(i).name;
    timestampStr = extractBetween(fileName, '_s', '_e');
    if isempty(timestampStr)
        continue;
    end
    timestampStr = timestampStr{1};
    timestampFormat = 'uuuuDDDHHmmssSSS'; 
    try
        fileTimestamp = datetime(timestampStr, 'InputFormat', timestampFormat);
        if fileTimestamp >= startDate && fileTimestamp <= endDate
            fileNames{end+1} = fileName;
            fileTimestamps(end+1) = fileTimestamp;
        end
    catch
        continue;
    end
end

if isempty(fileNames)
    error('No files found within the specified date range.');
end

[fileTimestamps, sortIdx] = sort(fileTimestamps);
fileNames = fileNames(sortIdx);

% Determine the variable name (CMI or Rad)
filePath = fullfile(downloadDir, fileNames{1});
info = ncinfo(filePath);
if any(strcmp({info.Variables.Name}, 'CMI'))
    variableName = 'CMI';
elseif any(strcmp({info.Variables.Name}, 'Rad'))
    variableName = 'Rad';
else
    availableVars = {info.Variables.Name};
    warning('No suitable variable found. Using first available variable: %s', availableVars{1});
    variableName = availableVars{1};
end

%% Process Frames (IR Data -> processedFrames)
processedFrames = cell(length(fileNames), 1);
for i = 1:length(fileNames)
    fileName = fileNames{i};
    filePath = fullfile(downloadDir, fileName);
    fprintf('Processing %s with method %s...\n', fileName, methodName);
    try
        data = double(ncread(filePath, variableName));
        
        % --- CHANGE #2: Replace NaNs with zero and print how many ---
        nanCount = sum(isnan(data(:)));
        if nanCount > 0
            fprintf('Replacing %d NaNs in %s\n', nanCount, fileName);
            data(isnan(data)) = 0;
        end
        
        data_preprocessed = preprocessData(data);
        img_processed = processDataMethod(data_preprocessed, methodName);
        processedFrames{i} = img_processed; % single-channel 2D array
    catch ME
        fprintf('Error processing %s: %s\n', fileName, ME.message);
        continue;
    end
end
fprintf('All frames processed using method: %s\n', methodName);

%% Parameters and Setup (Second Part)
% Instead of reading .nc files from VIS, we use processedFrames directly.

% Start timing the entire process (optional)
total_time = tic;

% Use all processed frames from the first part
num_frames = length(processedFrames);

if num_frames < 2
    error('At least two frames are required for processing.');
end

% Shrinking factor
shrinkfactor = 2;
invshrinkfactor = 1 / shrinkfactor;

% Spatial and wavelet parameters
degrees_per_pixel = 0.04;
km_per_degree = 111.32; 
original_pixel_size_km = degrees_per_pixel * km_per_degree;
pixel_size_km = original_pixel_size_km * shrinkfactor;

square_size_deg = 10; 
square_size_km = square_size_deg * km_per_degree;
brightness_threshold = 0.00001; 
std_threshold = 10;

% Wavelet parameters
NANGLES = 24;
Angles = 0:pi/(NANGLES-1):pi;
min_scale_km = 10;
max_scale_km = 500;
NSCALES=20;
Scales_km = logspace(log10(min_scale_km), log10(max_scale_km), NSCALES);
Scales = Scales_km / pixel_size_km;

window_buffer = 10; 
Preprocess_Flag = 1; 
time_interval = 1800; 
Save_Metadata_Flag = 0; 

% Define the parameters of the radial window
radius_factor = 1.1;  % 110% of the image radius (to keep the whole image)
decay_rate = 0.05;    % Control parameter for decay rate

%% Prepare First Frame for Dimension Calculation
% Take the first processed frame
frame1 = processedFrames{1};

% Resize if needed
if shrinkfactor ~= 1
    frame1 = imresize(frame1, invshrinkfactor);
end

% Preprocess first frame if needed
if Preprocess_Flag == 1
    frame1 = preprocess_img(frame1);
end

% Apply radial window after all other preprocessing is done
frame1_windowed = apply_radial_window(frame1, radius_factor, decay_rate);

% Now we have the final processed dimensions for the first frame
[frame_height_final, frame_width_final] = size(frame1_windowed);

% Adjust frame dimensions by removing buffer
x_buffer_range = (window_buffer + 1):(frame_width_final - window_buffer);
y_buffer_range = (window_buffer + 1):(frame_height_final - window_buffer);

adjusted_frame_width = length(x_buffer_range);
adjusted_frame_height = length(y_buffer_range);

% Calculate number of squares
square_size_px = round(square_size_km / pixel_size_km);
num_squares_x = ceil(adjusted_frame_width / square_size_px);
num_squares_y = ceil(adjusted_frame_height / square_size_px);

% Generate Squares Coordinates based on final processed frame dimensions
squares = [];
idx = 1;
for i = 1:num_squares_y
    for j = 1:num_squares_x
        x_start = floor((j - 1) * adjusted_frame_width / num_squares_x) + 1;
        y_start = floor((i - 1) * adjusted_frame_height / num_squares_y) + 1;
        x_end = floor(j * adjusted_frame_width / num_squares_x);
        y_end = floor(i * adjusted_frame_height / num_squares_y);

        if x_end > x_start && y_end > y_start
            squares(idx).x_range = x_buffer_range(x_start:x_end);
            squares(idx).y_range = y_buffer_range(y_start:y_end);
            squares(idx).index = idx;
            idx = idx + 1;
        end
    end
end
%%
% Display Frame 1 with Squares

figure;
% Preprocess frames if needed
if Preprocess_Flag==1
    frame1 = preprocess_img(frame1);
end


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

%% Create a video of all frames with the grid overlay
outputVideoFile = 'all_frames_with_grid.mj2'; % Nom du fichier vidéo

v = VideoWriter(outputVideoFile, 'Motion JPEG 2000');
v.FrameRate = 2;  % Fréquence d'images, ajustez selon vos besoins
v.LosslessCompression = true; % Compression sans perte
open(v);

figure('Position', [100, 100, 800, 600]); % Ajustez si nécessaire

for f_idx = 1:num_frames
    img = processedFrames{f_idx};

    % Resize si nécessaire
    if shrinkfactor ~= 1
        img = imresize(img, invshrinkfactor);
    end

    % Prétraitement si nécessaire
    if Preprocess_Flag
        img = preprocess_img(img);
    end

    % Appliquer la fenêtre radiale
    img_windowed = apply_radial_window(img, radius_factor, decay_rate);

    % Appliquer la zone tampon
    img_windowed = img_windowed(y_buffer_range, x_buffer_range);

    % Affichage
    imshow(img_windowed, [0, 1]);
    colormap('gray');
    hold on;

    % Dessin de la grille
    for s_idx = 1:length(squares)
        rectangle('Position', [squares(s_idx).x_range(1)-window_buffer, ...
                               squares(s_idx).y_range(1)-window_buffer, ...
                               length(squares(s_idx).x_range), ...
                               length(squares(s_idx).y_range)], ...
                  'EdgeColor', 'r', 'LineWidth', 1);
        text(double(squares(s_idx).x_range(1)-window_buffer), ...
             double(squares(s_idx).y_range(1)-window_buffer - 5), ...
             sprintf('%d', squares(s_idx).index), 'Color', 'yellow', 'FontSize', 8);
    end

    % Ajout du titre : numéro de frame et date exacte
    % On suppose que fileTimestamps contient le timing exact de chaque frame
    frame_timestamp = fileTimestamps(f_idx);
    title(sprintf('Frame %d - %s', f_idx, datestr(frame_timestamp, 'yyyy-mm-dd HH:MM:SS')));

    hold off;

    % Capturer l'image affichée pour la vidéo
    frame_img = getframe(gcf);
    writeVideo(v, frame_img);
end

close(v);
fprintf('Video created: %s\n', outputVideoFile);


%%
% Initialize Data Arrays
total_entries = (num_frames - 1) * length(squares) * NSCALES;
peak_data = NaN(total_entries, 8);
entry_counter = 0;
square_metadata = struct('frame_idx', {}, 'square_idx', {}, 'spec1', {});

% Process Frames in Pairs
for frame_idx = 1:num_frames - 1
    fprintf('Processing frame pair %d/%d...\n', frame_idx, num_frames);
    frame_time = tic;

    % Get current and next frame from processedFrames
    frame1 = processedFrames{frame_idx};
    frame2 = processedFrames{frame_idx + 1};

    % Transpose if needed
    %frame1 = frame1';
    %frame2 = frame2';

    % Resize if needed
    if shrinkfactor ~= 1
        frame1 = imresize(frame1, invshrinkfactor);
        frame2 = imresize(frame2, invshrinkfactor);
    end

    % Preprocess frames if needed
    if Preprocess_Flag
        frame1 = preprocess_img(frame1);
        frame2 = preprocess_img(frame2);
    end

    % Apply radial window
    frame1_windowed = apply_radial_window(frame1, radius_factor, decay_rate);
    frame2_windowed = apply_radial_window(frame2, radius_factor, decay_rate);

    % Apply the previously defined buffer ranges
    frame1_windowed = frame1_windowed(y_buffer_range, x_buffer_range);
    frame2_windowed = frame2_windowed(y_buffer_range, x_buffer_range);

    % Compute CWT for both frames
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

        if Save_Metadata_Flag == 1
            metadata_entry.frame_idx = frame_idx;
            metadata_entry.square_idx = s_idx;
            metadata_entry.spec1 = spec1;
            square_metadata = [square_metadata; metadata_entry];
        end

        img1 = frame1_windowed(y_range, x_range);
        mean_brightness = mean(img1(:));
        std_brightness = std(img1(:));

        if mean_brightness < brightness_threshold || std_brightness > std_threshold
            continue; % Skip this square if it doesn't meet criteria
        end

        % Compute cross-wavelet spectrum
        xwt = spec1 .* conj(spec2);
        coherence = abs(xwt);
        phase_difference = angle(xwt);

        % Find peaks and speeds
        peak_list = find_peaks_and_speeds(coherence, phase_difference, Scales, Angles, pixel_size_km, time_interval);
        n_scales = size(peak_list, 1);

        data_to_store = [repmat(frame_idx, n_scales, 1), repmat(s_idx, n_scales, 1), peak_list];
        idx_start = entry_counter + 1;
        idx_end = entry_counter + n_scales;

        if idx_end > total_entries
            extra_entries = total_entries; % double size if needed
            peak_data = [peak_data; NaN(extra_entries, 8)];
            total_entries = total_entries + extra_entries;
        end

        peak_data(idx_start:idx_end, :) = data_to_store;
        entry_counter = idx_end;
    end

    fprintf('Time to process frame pair %d: %.2f seconds.\n', frame_idx, toc(frame_time));
end

peak_data = peak_data(1:entry_counter, :);
fprintf('Total execution time: %.2f seconds.\n', toc(total_time));
%% Analyze the data

%plot_waverose(1,1,false,4);

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
    
    % Truncate values below the 10th percentile to the lower bound
    % and values above the 99th percentile to the upper bound.
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
        %wavelength_km = Scales(scale_idx) * 2 * pixel_size_km;  % Wavelength in km
        wavelength_km = Scales(scale_idx) *pi/sqrt(2) * pixel_size_km;  % Wavelength in km
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

function plot_waverose_old(frame_id, square_id, display_sine_wave_plots, figs_per_page)
%% Plot Waverose Function with Enhanced Features
% This function generates an advanced rose plot for a specific frame and square.
% It uses data already processed and stored in the base workspace.
%
% Parameters:
% - frame_id: ID of the frame to process (index in processedFrames)
% - square_id: ID of the square within the frame (from squares struct)
% - display_sine_wave_plots: (optional) true/false to display sine wave plots (default: true)
% - figs_per_page: (optional) number of subplots per figure page (default: 20)

% Set default values for optional parameters
if nargin < 3
    display_sine_wave_plots = true;
end
if nargin < 4
    figs_per_page = 20;
end

% Retrieve necessary variables from the base workspace
processedFrames     = evalin('base', 'processedFrames');
invshrinkfactor     = evalin('base', 'invshrinkfactor');
pixel_size_km       = evalin('base', 'pixel_size_km');
Scales              = evalin('base', 'Scales');
Angles              = evalin('base', 'Angles');
time_interval       = evalin('base', 'time_interval');
squares             = evalin('base', 'squares');
Preprocess_Flag     = evalin('base', 'Preprocess_Flag');
window_buffer       = evalin('base', 'window_buffer');
radius_factor       = evalin('base', 'radius_factor');
decay_rate          = evalin('base', 'decay_rate');
shrinkfactor        = evalin('base', 'shrinkfactor');
brightness_threshold= evalin('base', 'brightness_threshold');
std_threshold       = evalin('base', 'std_threshold');

% Validate frame_id
num_frames = length(processedFrames);
if frame_id < 1 || frame_id >= num_frames
    error('Invalid frame_id. Must be between 1 and %d.', num_frames - 1);
end

% Retrieve frames directly from processedFrames
frame1 = processedFrames{frame_id};
frame2 = processedFrames{frame_id + 1};

% Resize if needed
if shrinkfactor ~= 1
    frame1 = imresize(frame1, invshrinkfactor);
    frame2 = imresize(frame2, invshrinkfactor);
end

% Preprocess frames if needed
if Preprocess_Flag
    frame1 = preprocess_img(frame1);
    frame2 = preprocess_img(frame2);
end

% Apply radial window
frame1_windowed = apply_radial_window(frame1, radius_factor, decay_rate);
frame2_windowed = apply_radial_window(frame2, radius_factor, decay_rate);

% Retrieve buffer ranges (these should already be in the base workspace)
x_buffer_range = evalin('base', 'x_buffer_range');
y_buffer_range = evalin('base', 'y_buffer_range');

% Apply buffer ranges
frame1_windowed = frame1_windowed(y_buffer_range, x_buffer_range);
frame2_windowed = frame2_windowed(y_buffer_range, x_buffer_range);

% Perform Wavelet Transform on Full Frames
cwt1_full = cwtft2(frame1_windowed, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
spec1_full = squeeze(cwt1_full.cfs);

cwt2_full = cwtft2(frame2_windowed, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
spec2_full = squeeze(cwt2_full.cfs);

%% Extract the Specified Square
% Find the square with the given square_id
idx = find([squares.index] == square_id);
if isempty(idx)
    error('Invalid square_id. Square not found.');
end
square = squares(idx);

% Use the square's x_range and y_range directly
x_range = square.x_range;
y_range = square.y_range;

% Ensure x_range, y_range fit in spec arrays
x_range(x_range < 1) = [];
y_range(y_range < 1) = [];
x_range(x_range > size(spec1_full, 2)) = [];
y_range(y_range > size(spec1_full, 1)) = [];

% Extract the cwt coefficients for the square
spec1 = spec1_full(y_range, x_range, :, :);
spec2 = spec2_full(y_range, x_range, :, :);

% Compute power spectra of spec1 and spec2
power1 = abs(spec1).^2;

% Check brightness and standard deviation on the original frames (frame1_windowed)
img1 = frame1_windowed(y_range, x_range);
mean_brightness = mean(img1(:));
std_brightness = std(img1(:));

if mean_brightness < brightness_threshold || std_brightness > std_threshold
    error('The selected square does not meet brightness or standard deviation thresholds.');
end

%% Compute Cross-Wavelet Spectrum, Coherence, Phase Difference
xwt = spec1 .* conj(spec2);
coherence = abs(xwt);
phase_difference = angle(xwt);

%% Find Peaks and Speeds
peak_list = find_peaks_and_speeds(coherence, phase_difference, Scales, Angles, pixel_size_km, time_interval);

%% Plot Advanced Rose Plot
% Generate the advanced rose plot with power and coherence
% Overlay the peaks on the plot

% Calculate inner power and coherence for plotting
buffer = 0; % Adjust buffer as needed
innerpower = squeeze(mean(mean(power1(buffer+1:end-buffer, buffer+1:end-buffer, :, :), 1, 'omitnan'), 2, 'omitnan')); % (num_scales x num_angles)
innercoherence = squeeze(mean(mean(coherence(buffer+1:end-buffer, buffer+1:end-buffer, :, :), 1, 'omitnan'), 2, 'omitnan')); % (num_scales x num_angles)

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
colormap(ax2, 'autumn');
axis equal;
set(ax2, 'Position', [0.1, 0.1, 0.75, 0.75]);
ax2.XTick = [];
ax2.YTick = [];
set(ax2, 'Color', 'none');
linkaxes([ax1, ax2]);
hold on;

% Overlay peaks on the coherence plot
if ~isempty(peak_list)
    % Since peak_list is now (num_scales x 8)
    max_scales = peak_list(:,1); % Scales
    max_angles = peak_list(:,2) + pi;  % Adjust angles for lower half
    [peak_X, peak_Y] = pol2cart(max_angles, max_scales);
    plot(ax2, peak_X, peak_Y, 'k*', 'MarkerSize', 10);
    % Annotate peaks with speeds 
    for i = 1:length(peak_X)
        % Display speed in m/s
        text(ax2, peak_X(i) * 1.05, peak_Y(i) * 1.05, sprintf('%.1f m/s', peak_list(i, 4)), 'Color', 'k', 'FontSize', 10);
    end
end

% Adjust axes limits
xlim(ax1, [min(X_pos(:)) - 1, max(X_pos(:)) + 1]);
ylim(ax1, [min(Y_neg(:)) - 1, max(Y_pos(:)) + 1]);

%% Add Radial Rings and Angle Labels
% Add radial rings corresponding to scales (in km)
wavelengths_km = Scales * pi/sqrt(2) * pixel_size_km; % Convert scales to wavelengths in km
ring_radii = Scales;
for i = 1:length(ring_radii)
    theta_ring = linspace(0, 2 * pi, 100);
    [x_ring, y_ring] = pol2cart(theta_ring, ring_radii(i));
    plot(ax1, x_ring, y_ring, 'k--');
    plot(ax2, x_ring, y_ring, 'k--');
    % Add scale labels in km
    text(ax1, ring_radii(i) * 1.05, 0, sprintf('%.1f km', wavelengths_km(i)), 'HorizontalAlignment', 'left');
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

%% Additional Plots: Coherence Maps and Phase Difference Maps

% Define common fractions of pi for displaying angles
pi_fractions = {'0', '\pi/6', '\pi/4', '\pi/3', '\pi/2', '2\pi/3', '\pi', '4\pi/3', '3\pi/2', '2\pi'};
pi_fraction_values = [0, pi/6, pi/4, pi/3, pi/2, 2*pi/3, pi, 4*pi/3, 3*pi/2, 2*pi];

% Preallocate a cell array to store coherence masks
coherence_masks = cell(size(peak_list,1),1);

% Plot coherence with mask for each scale
num_peaks = size(peak_list, 1);

% Calculate the number of figures needed based on figs_per_page
num_figures = ceil(num_peaks / figs_per_page);
current_peak = 1;

for fig_num = 1:num_figures
    figure;
    num_subplots = min(figs_per_page, num_peaks - (fig_num - 1) * figs_per_page);
    ncols = ceil(sqrt(num_subplots));
    nrows = ceil(num_subplots / ncols);
    
    for subplot_idx = 1:num_subplots
        i = current_peak;
        
        % Extract the real scale and angle from the peak_list
        real_scale = peak_list(i,1);  % Scale from peak_list
        real_angle = peak_list(i,2);  % Angle from peak_list
        
        % Find the index of the closest scale in the Scales array
        [~, scale_idx] = min(abs(Scales - real_scale));
        
        % Adjust real_angle if necessary (wrap around)
        real_angle = mod(real_angle, 2*pi);
        
        % Find the index of the closest angle in the Angles array
        [~, angle_idx] = min(abs(Angles - real_angle));
        
        % Extract the corresponding coherence slice for the current scale and angle
        coherence_slice = coherence(:, :, scale_idx, angle_idx);
        
        % Define the threshold for the top 60% of the max coherence value
        coherence_max = max(coherence_slice(:));
        coherence_threshold = 0.6 * coherence_max;
        
        % Create the mask where coherence is above the threshold
        coherence_mask = coherence_slice >= coherence_threshold;
        
        % Store the mask in the cell array
        coherence_masks{i} = coherence_mask;  % Store mask for later use
        
        % Plot the coherence using imagesc
        subplot(nrows, ncols, subplot_idx);  % Subplot for multiple plots
        imagesc(coherence_slice);
        hold on;
        % Overlay the contour of the coherence mask
        contour(coherence_mask, [1 1], 'r', 'LineWidth', 1);  % Red contour at mask boundary
        
        % Ensure axis is tight
        axis tight;
        
        % Convert the real_angle to a fraction of pi for the title
        [~, angle_fraction_idx] = min(abs(pi_fraction_values - real_angle));  % Find closest pi fraction
        angle_str = pi_fractions{angle_fraction_idx};  % Get the corresponding fraction of pi string
        
        % Set title with scale and angle in fractions of pi
        wavelength_km = real_scale * pi/sqrt(2) * pixel_size_km;
        title(sprintf('Wavelength: %.1f km, Angle: %s', wavelength_km, angle_str));
        
        % Customize the colorbar
        colorbar;
        
        % Set axis equal for consistent plotting
        axis equal;
        hold off;
        
        current_peak = current_peak + 1;
    end
    
    sgtitle('Coherence with Mask for Each Scale/Angle Peak');
end

%% Plot Phase Difference with Masked Values and Compute Mean
current_peak = 1;
for fig_num = 1:num_figures
    figure;
    num_subplots = min(figs_per_page, num_peaks - (fig_num - 1) * figs_per_page);
    ncols = ceil(sqrt(num_subplots));
    nrows = ceil(num_subplots / ncols);
    
    for subplot_idx = 1:num_subplots
        i = current_peak;
        
        % Extract the real scale and angle from the peak_list
        real_scale = peak_list(i,1);  % Scale from peak_list
        real_angle = peak_list(i,2);  % Angle from peak_list
        
        % Find the index of the closest scale in the Scales array
        [~, scale_idx] = min(abs(Scales - real_scale));
        
        % Adjust real_angle if necessary (wrap around)
        real_angle = mod(real_angle, 2*pi);
        
        % Find the index of the closest angle in the Angles array
        [~, angle_idx] = min(abs(Angles - real_angle));
        
        % Convert the real_angle to a fraction of pi for the title
        [~, angle_fraction_idx] = min(abs(pi_fraction_values - real_angle));  % Find closest pi fraction
        angle_str = pi_fractions{angle_fraction_idx};  % Get the corresponding fraction of pi string
        
        % Extract the corresponding phase_difference slice for the current scale and angle
        phase_slice = phase_difference(:, :, scale_idx, angle_idx);
        
        % Retrieve the corresponding mask from coherence_masks
        coherence_mask = coherence_masks{i};
        
        % Apply the mask from coherence to the phase slice
        phase_masked = phase_slice;
        phase_masked(~coherence_mask) = NaN;  % Set non-coherent areas to NaN
        
        % Plot the phase_difference using imagesc
        subplot(nrows, ncols, subplot_idx);  % Subplot for multiple plots
        imagesc(phase_masked);
        hold on;
        
        % Overlay the contour of the coherence mask
        contour(coherence_mask, [1 1], 'r', 'LineWidth', 1);  % Red contour at mask boundary
        
        % Set title with scale and angle in fractions of pi
        wavelength_km = real_scale * pi/sqrt(2) * pixel_size_km;
        title(sprintf('Wavelength: %.1f km, Angle: %s', wavelength_km, angle_str));
        
        % Customize the colorbar to display ticks from -pi to pi
        c = colorbar;
        caxis([-pi pi]);  % Set color axis limits from -pi to pi
        set(c, 'Ticks', [-pi, -pi/2, 0, pi/2, pi], 'TickLabels', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
        
        % Set axis equal for consistent plotting
        axis equal;
        axis tight;
        
        hold off;
        
        current_peak = current_peak + 1;
    end
    
    sgtitle('Phase Difference with Masked Regions for Each Scale/Angle Peak');
end

%% Additional Plot: Overlay Wavelet Contours on the Image Zoomed into the Square
current_peak = 1;
for fig_num = 1:num_figures
    figure;
    num_subplots = min(figs_per_page, num_peaks - (fig_num - 1) * figs_per_page);
    ncols = ceil(sqrt(num_subplots));
    nrows = ceil(num_subplots / ncols);
    
    for subplot_idx = 1:num_subplots
        i = current_peak;
        
        % Extract the real scale and angle from the peak_list
        real_scale = peak_list(i,1);  % Scale from peak_list
        real_angle = peak_list(i,2);  % Angle from peak_list
        
        % Find the index of the closest scale in the Scales array
        [~, scale_idx] = min(abs(Scales - real_scale));
        
        % Adjust real_angle if necessary (wrap around)
        real_angle = mod(real_angle, 2*pi);
        
        % Find the index of the closest angle in the Angles array
        [~, angle_idx] = min(abs(Angles - real_angle));
        
        % Extract the corresponding frame square (use the windowed frame)
        frame_square = frame1_windowed(y_range, x_range);
        
        % Extract the wavelet coefficients for the square
        spec_square = spec1(:, :, :, :);
        
        % Plot using image_with_wavelet_overlay
        subplot(nrows, ncols, subplot_idx);
        clevfactor = 1;  % Adjust as needed
        ProcessFlag = 1; % Use ProcessFlag = 1 for normalized display
        
        % Call the overlay function
        image_with_wavelet_overlay(frame_square, spec_square, Scales, scale_idx, angle_idx, clevfactor, ProcessFlag);
        
        % Overlay the mask contour from coherence_mask
        coherence_mask = coherence_masks{i};  % Get the mask for this peak
        hold on;
        contour(coherence_mask, [1 1], 'magenta', 'LineWidth', 1);  % Magenta contour at mask boundary
        hold off;
        
        % Set title with scale and angle
        wavelength_km = real_scale * pi/sqrt(2) * pixel_size_km;
        title(sprintf('Wavelength: %.1f km, Angle: %.2f°', wavelength_km, real_angle * (180/pi)));
        
        % Ensure axis is equal and tight
        axis equal;
        axis tight;
        
        current_peak = current_peak + 1;
    end
    
    sgtitle('Wavelet Contours Overlaid on Image Zoomed into Square');
end

%% Calculate Speed from Phase Shift

% Initialize arrays to store results
speeds = zeros(num_peaks, 1);

for i = 1:num_peaks
    % Extract scale and mean phase difference from peak_list
    scale = peak_list(i,1);  % Scale in pixels
    mean_phase_difference = peak_list(i,3);  % Mean phase difference in radians
    
    % Calculate the wavelength in km
    wavelength_km = scale * pi/sqrt(2) * pixel_size_km;
    
    % Calculate the distance shift in km
    distance_shift_km = mean_phase_difference * wavelength_km / (2 * pi);
    
    % Calculate the speed in km/s
    speed_km_per_s = distance_shift_km / time_interval;
    
    % Store the results
    speeds(i) = speed_km_per_s;
    
    % Update the fourth column of peak_list with the speed in m/s
    peak_list(i,4) = speed_km_per_s * 1000;  % Convert to m/s
    
    % Only display sine wave plots if the option is true
    if display_sine_wave_plots
        % Plot the sine waves (optional)
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
        title(sprintf('Sine Wave with Phase Shift for Peak %d\nWavelength: %.2f km\nPhase: %.4f radians (%.4f\\pi)\nShift: %.4f km, Speed: %.4f m/s', ...
            i, wavelength_km, mean_phase_difference, mean_phase_difference / pi, distance_shift_km, speed_km_per_s * 1000));
        grid on;
        xlim([0, wavelength_km]);
        hold off;
    end
end

% Display the updated peak_list
%disp('Updated peak_list with mean phase values and speeds:');
%disp('Columns: 1-Scale, 2-Angle, 3-Mean Phase Difference, 4-Speed (m/s), 5-Coherence Value, 6-Wavelength_km');
%disp(peak_list);

end

function plot_waverose(frame_id, square_id, display_sine_wave_plots, figs_per_page)
%% Plot Waverose Function with Enhanced Features
% This function generates an advanced rose plot for a specific frame and square.
% It uses data already processed and stored in the base workspace.
%
% Parameters:
% - frame_id: ID of the frame to process (index in processedFrames)
% - square_id: ID of the square within the frame (from squares struct)
% - display_sine_wave_plots: (optional) true/false to display sine wave plots (default: true)
% - figs_per_page: (optional) number of subplots per figure page (default: 20)

% Set default values for optional parameters
if nargin < 3
    display_sine_wave_plots = true;
end
if nargin < 4
    figs_per_page = 20;
end

% Retrieve necessary variables from the base workspace
processedFrames     = evalin('base', 'processedFrames');
invshrinkfactor     = evalin('base', 'invshrinkfactor');
pixel_size_km       = evalin('base', 'pixel_size_km');
Scales              = evalin('base', 'Scales');
Angles              = evalin('base', 'Angles');
time_interval       = evalin('base', 'time_interval');
squares             = evalin('base', 'squares');
Preprocess_Flag     = evalin('base', 'Preprocess_Flag');
window_buffer       = evalin('base', 'window_buffer');
radius_factor       = evalin('base', 'radius_factor');
decay_rate          = evalin('base', 'decay_rate');
shrinkfactor        = evalin('base', 'shrinkfactor');
brightness_threshold= evalin('base', 'brightness_threshold');
std_threshold       = evalin('base', 'std_threshold');

% Validate frame_id
num_frames = length(processedFrames);
if frame_id < 1 || frame_id >= num_frames
    error('Invalid frame_id. Must be between 1 and %d.', num_frames - 1);
end

% Retrieve frames directly from processedFrames
frame1 = processedFrames{frame_id};
frame2 = processedFrames{frame_id + 1};

% Resize if needed
if shrinkfactor ~= 1
    frame1 = imresize(frame1, invshrinkfactor);
    frame2 = imresize(frame2, invshrinkfactor);
end

% Preprocess frames if needed
if Preprocess_Flag
    frame1 = preprocess_img(frame1);
    frame2 = preprocess_img(frame2);
end

% Apply radial window
frame1_windowed = apply_radial_window(frame1, radius_factor, decay_rate);
frame2_windowed = apply_radial_window(frame2, radius_factor, decay_rate);

% Retrieve buffer ranges (these should already be in the base workspace)
x_buffer_range = evalin('base', 'x_buffer_range');
y_buffer_range = evalin('base', 'y_buffer_range');

% Apply buffer ranges
frame1_windowed = frame1_windowed(y_buffer_range, x_buffer_range);
frame2_windowed = frame2_windowed(y_buffer_range, x_buffer_range);

% Perform Wavelet Transform on Full Frames
cwt1_full = cwtft2(frame1_windowed, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
spec1_full = squeeze(cwt1_full.cfs);

cwt2_full = cwtft2(frame2_windowed, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
spec2_full = squeeze(cwt2_full.cfs);

%% Extract the Specified Square
% Find the square with the given square_id
idx = find([squares.index] == square_id);
if isempty(idx)
    error('Invalid square_id. Square not found.');
end
square = squares(idx);

% Use the square's x_range and y_range directly
x_range = square.x_range;
y_range = square.y_range;

% Ensure x_range, y_range fit in spec arrays
x_range(x_range < 1) = [];
y_range(y_range < 1) = [];
x_range(x_range > size(spec1_full, 2)) = [];
y_range(y_range > size(spec1_full, 1)) = [];

% Extract the cwt coefficients for the square
spec1 = spec1_full(y_range, x_range, :, :);
spec2 = spec2_full(y_range, x_range, :, :);

% Check brightness and standard deviation on the original frames (frame1_windowed)
img1 = frame1_windowed(y_range, x_range);
mean_brightness = mean(img1(:));
std_brightness = std(img1(:));

if mean_brightness < brightness_threshold || std_brightness > std_threshold
    error('The selected square does not meet brightness or standard deviation thresholds.');
end

%% Compute Cross-Wavelet Spectrum, Coherence, Phase Difference
xwt = spec1 .* conj(spec2);
coherence = abs(xwt);
phase_difference = angle(xwt);

for s_idx = 1:length(Scales)
    this_abs_xwt = coherence(:,:,s_idx,:);
    scale_mean = mean(this_abs_xwt(:), 'omitnan'); 
    if scale_mean ~= 0
        coherence(:,:,s_idx,:) = coherence(:,:,s_idx,:) / scale_mean;
    end
end


% Compute power spectra of spec1 and spec2
power1 = abs(spec1).^2;

for s_idx = 1:length(Scales)
    this_abs_xwt = power1(:,:,s_idx,:);
    scale_mean = sqrt(mean(this_abs_xwt(:).^2, 'omitnan')); 
    if scale_mean ~= 0
        power1(:,:,s_idx,:) = power1(:,:,s_idx,:) / scale_mean;
    end
end


%% Find Peaks and Speeds
peak_list = find_peaks_and_speeds(coherence, phase_difference, Scales, Angles, pixel_size_km, time_interval);

%% Plot Advanced Rose Plot
% Generate the advanced rose plot with power and coherence
% Overlay the peaks on the plot

% Calculate inner power and coherence for plotting
buffer = 0; % Adjust buffer as needed
innerpower = squeeze(mean(mean(power1(buffer+1:end-buffer, buffer+1:end-buffer, :, :), 1, 'omitnan'), 2, 'omitnan')); % (num_scales x num_angles)
innercoherence = squeeze(mean(mean(coherence(buffer+1:end-buffer, buffer+1:end-buffer, :, :), 1, 'omitnan'), 2, 'omitnan')); % (num_scales x num_angles)

% Normalize power and coherence by scale (already done in the previous step)
% mean_power_byscale = mean(innerpower, 2);
% mean_coherence_byscale = mean(innercoherence, 2);
% 
% anglespec_power = innerpower ./ mean_power_byscale;
% anglespec_coherence = innercoherence ./ mean_coherence_byscale;

anglespec_power = innerpower ;
anglespec_coherence = innercoherence;

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
colormap(ax2, 'autumn');
axis equal;
set(ax2, 'Position', [0.1, 0.1, 0.75, 0.75]);
ax2.XTick = [];
ax2.YTick = [];
set(ax2, 'Color', 'none');
linkaxes([ax1, ax2]);
hold on;

% Overlay peaks on the coherence plot
if ~isempty(peak_list)
    % Since peak_list is now (num_scales x 8)
    max_scales = peak_list(:,1); % Scales
    max_angles = peak_list(:,2) + pi;  % Adjust angles for lower half
    [peak_X, peak_Y] = pol2cart(max_angles, max_scales);
    plot(ax2, peak_X, peak_Y, 'k*', 'MarkerSize', 10);
    % Annotate peaks with speeds 
    for i = 1:length(peak_X)
        % Display speed in m/s
        text(ax2, peak_X(i) * 1.05, peak_Y(i) * 1.05, sprintf('%.1f m/s', peak_list(i, 4)), 'Color', 'k', 'FontSize', 10);
    end
end

% Adjust axes limits
xlim(ax1, [min(X_pos(:)) - 1, max(X_pos(:)) + 1]);
ylim(ax1, [min(Y_neg(:)) - 1, max(Y_pos(:)) + 1]);

%% Add Radial Rings and Angle Labels
% Add radial rings corresponding to scales (in km)
wavelengths_km = Scales * pi/sqrt(2) * pixel_size_km; % Convert scales to wavelengths in km
ring_radii = Scales;
for i = 1:length(ring_radii)
    theta_ring = linspace(0, 2 * pi, 100);
    [x_ring, y_ring] = pol2cart(theta_ring, ring_radii(i));
    plot(ax1, x_ring, y_ring, 'k--');
    plot(ax2, x_ring, y_ring, 'k--');
    % Add scale labels in km
    text(ax1, ring_radii(i) * 1.05, 0, sprintf('%.1f km', wavelengths_km(i)), 'HorizontalAlignment', 'left');
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

%% Additional Plots: Coherence Maps and Phase Difference Maps

% Define common fractions of pi for displaying angles
pi_fractions = {'0', '\pi/6', '\pi/4', '\pi/3', '\pi/2', '2\pi/3', '\pi', '4\pi/3', '3\pi/2', '2\pi'};
pi_fraction_values = [0, pi/6, pi/4, pi/3, pi/2, 2*pi/3, pi, 4*pi/3, 3*pi/2, 2*pi];

% Preallocate a cell array to store coherence masks
coherence_masks = cell(size(peak_list,1),1);

% Plot coherence with mask for each scale
num_peaks = size(peak_list, 1);

% Calculate the number of figures needed based on figs_per_page
num_figures = ceil(num_peaks / figs_per_page);
current_peak = 1;

for fig_num = 1:num_figures
    figure;
    num_subplots = min(figs_per_page, num_peaks - (fig_num - 1) * figs_per_page);
    ncols = ceil(sqrt(num_subplots));
    nrows = ceil(num_subplots / ncols);
    
    for subplot_idx = 1:num_subplots
        i = current_peak;
        
        % Extract the real scale and angle from the peak_list
        real_scale = peak_list(i,1);  % Scale from peak_list
        real_angle = peak_list(i,2);  % Angle from peak_list
        
        % Find the index of the closest scale in the Scales array
        [~, scale_idx] = min(abs(Scales - real_scale));
        
        % Adjust real_angle if necessary (wrap around)
        real_angle = mod(real_angle, 2*pi);
        
        % Find the index of the closest angle in the Angles array
        [~, angle_idx] = min(abs(Angles - real_angle));
        
        % Extract the corresponding coherence slice for the current scale and angle
        coherence_slice = coherence(:, :, scale_idx, angle_idx);
        
        % Define the threshold for the top 60% of the max coherence value
        coherence_max = max(coherence_slice(:));
        coherence_threshold = 0.6 * coherence_max;
        
        % Create the mask where coherence is above the threshold
        coherence_mask = coherence_slice >= coherence_threshold;
        
        % Store the mask in the cell array
        coherence_masks{i} = coherence_mask;  % Store mask for later use
        
        % Plot the coherence using imagesc
        subplot(nrows, ncols, subplot_idx);  % Subplot for multiple plots
        imagesc(coherence_slice);
        hold on;
        % Overlay the contour of the coherence mask
        contour(coherence_mask, [1 1], 'r', 'LineWidth', 1);  % Red contour at mask boundary
        
        % Ensure axis is tight
        axis tight;
        
        % Convert the real_angle to a fraction of pi for the title
        [~, angle_fraction_idx] = min(abs(pi_fraction_values - real_angle));  % Find closest pi fraction
        angle_str = pi_fractions{angle_fraction_idx};  % Get the corresponding fraction of pi string
        
        % Set title with scale and angle in fractions of pi
        wavelength_km = real_scale * pi/sqrt(2) * pixel_size_km;
        title(sprintf('Wavelength: %.1f km, Angle: %s', wavelength_km, angle_str));
        
        % Customize the colorbar
        colorbar;
        
        % Set axis equal for consistent plotting
        axis equal;
        hold off;
        
        current_peak = current_peak + 1;
    end
    
    sgtitle('Coherence with Mask for Each Scale/Angle Peak');
end

%% Plot Phase Difference with Masked Values and Compute Mean
current_peak = 1;
for fig_num = 1:num_figures
    figure;
    num_subplots = min(figs_per_page, num_peaks - (fig_num - 1) * figs_per_page);
    ncols = ceil(sqrt(num_subplots));
    nrows = ceil(num_subplots / ncols);
    
    for subplot_idx = 1:num_subplots
        i = current_peak;
        
        % Extract the real scale and angle from the peak_list
        real_scale = peak_list(i,1);  % Scale from peak_list
        real_angle = peak_list(i,2);  % Angle from peak_list
        
        % Find the index of the closest scale in the Scales array
        [~, scale_idx] = min(abs(Scales - real_scale));
        
        % Adjust real_angle if necessary (wrap around)
        real_angle = mod(real_angle, 2*pi);
        
        % Find the index of the closest angle in the Angles array
        [~, angle_idx] = min(abs(Angles - real_angle));
        
        % Convert the real_angle to a fraction of pi for the title
        [~, angle_fraction_idx] = min(abs(pi_fraction_values - real_angle));  % Find closest pi fraction
        angle_str = pi_fractions{angle_fraction_idx};  % Get the corresponding fraction of pi string
        
        % Extract the corresponding phase_difference slice for the current scale and angle
        phase_slice = phase_difference(:, :, scale_idx, angle_idx);
        
        % Retrieve the corresponding mask from coherence_masks
        coherence_mask = coherence_masks{i};
        
        % Apply the mask from coherence to the phase slice
        phase_masked = phase_slice;
        phase_masked(~coherence_mask) = NaN;  % Set non-coherent areas to NaN
        
        % Plot the phase_difference using imagesc
        subplot(nrows, ncols, subplot_idx);  % Subplot for multiple plots
        imagesc(phase_masked);
        hold on;
        
        % Overlay the contour of the coherence mask
        contour(coherence_mask, [1 1], 'r', 'LineWidth', 1);  % Red contour at mask boundary
        
        % Set title with scale and angle in fractions of pi
        wavelength_km = real_scale * pi/sqrt(2) * pixel_size_km;
        title(sprintf('Wavelength: %.1f km, Angle: %s', wavelength_km, angle_str));
        
        % Customize the colorbar to display ticks from -pi to pi
        c = colorbar;
        caxis([-pi pi]);  % Set color axis limits from -pi to pi
        set(c, 'Ticks', [-pi, -pi/2, 0, pi/2, pi], 'TickLabels', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
        
        % Set axis equal for consistent plotting
        axis equal;
        axis tight;
        
        hold off;
        
        current_peak = current_peak + 1;
    end
    
    sgtitle('Phase Difference with Masked Regions for Each Scale/Angle Peak');
end

%% Additional Plot: Overlay Wavelet Contours on the Image Zoomed into the Square
current_peak = 1;
for fig_num = 1:num_figures
    figure;
    num_subplots = min(figs_per_page, num_peaks - (fig_num - 1) * figs_per_page);
    ncols = ceil(sqrt(num_subplots));
    nrows = ceil(num_subplots / ncols);
    
    for subplot_idx = 1:num_subplots
        i = current_peak;
        
        % Extract the real scale and angle from the peak_list
        real_scale = peak_list(i,1);  % Scale from peak_list
        real_angle = peak_list(i,2);  % Angle from peak_list
        
        % Find the index of the closest scale in the Scales array
        [~, scale_idx] = min(abs(Scales - real_scale));
        
        % Adjust real_angle if necessary (wrap around)
        real_angle = mod(real_angle, 2*pi);
        
        % Find the index of the closest angle in the Angles array
        [~, angle_idx] = min(abs(Angles - real_angle));
        
        % Extract the corresponding frame square (use the windowed frame)
        frame_square = frame1_windowed(y_range, x_range);
        
        % Extract the wavelet coefficients for the square
        spec_square = spec1(:, :, :, :);
        
        % Plot using image_with_wavelet_overlay
        subplot(nrows, ncols, subplot_idx);
        clevfactor = 1;  % Adjust as needed
        ProcessFlag = 1; % Use ProcessFlag = 1 for normalized display
        
        % Call the overlay function
        image_with_wavelet_overlay(frame_square, spec_square, Scales, scale_idx, angle_idx, clevfactor, ProcessFlag);
        
        % Overlay the mask contour from coherence_mask
        coherence_mask = coherence_masks{i};  % Get the mask for this peak
        hold on;
        contour(coherence_mask, [1 1], 'magenta', 'LineWidth', 1);  % Magenta contour at mask boundary
        hold off;
        
        % Set title with scale and angle
        wavelength_km = real_scale * pi/sqrt(2) * pixel_size_km;
        title(sprintf('Wavelength: %.1f km, Angle: %.2f°', wavelength_km, real_angle * (180/pi)));
        
        % Ensure axis is equal and tight
        axis equal;
        axis tight;
        
        current_peak = current_peak + 1;
    end
    
    sgtitle('Wavelet Contours Overlaid on Image Zoomed into Square');
end

%% Calculate Speed from Phase Shift

% Initialize arrays to store results
speeds = zeros(num_peaks, 1);

for i = 1:num_peaks
    % Extract scale and mean phase difference from peak_list
    scale = peak_list(i,1);  % Scale in pixels
    mean_phase_difference = peak_list(i,3);  % Mean phase difference in radians
    
    % Calculate the wavelength in km
    wavelength_km = scale * pi/sqrt(2) * pixel_size_km;
    
    % Calculate the distance shift in km
    distance_shift_km = mean_phase_difference * wavelength_km / (2 * pi);
    
    % Calculate the speed in km/s
    speed_km_per_s = distance_shift_km / time_interval;
    
    % Store the results
    speeds(i) = speed_km_per_s;
    
    % Update the fourth column of peak_list with the speed in m/s
    peak_list(i,4) = speed_km_per_s * 1000;  % Convert to m/s
    
    % Only display sine wave plots if the option is true
    if display_sine_wave_plots
        % Plot the sine waves (optional)
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
        title(sprintf('Sine Wave with Phase Shift for Peak %d\nWavelength: %.2f km\nPhase: %.4f radians (%.4f\\pi)\nShift: %.4f km, Speed: %.4f m/s', ...
            i, wavelength_km, mean_phase_difference, mean_phase_difference / pi, distance_shift_km, speed_km_per_s * 1000));
        grid on;
        xlim([0, wavelength_km]);
        hold off;
    end
end

% Display the updated peak_list
%disp('Updated peak_list with mean phase values and speeds:');
%disp('Columns: 1-Scale, 2-Angle, 3-Mean Phase Difference, 4-Speed (m/s), 5-Coherence Value, 6-Wavelength_km');
%disp(peak_list);

end

function image_with_wavelet_overlay(img, spec, Scales, scale_idx, angle_idx, clevfactor, ProcessFlag)
    % Function to overlay wavelet contours on an image
    %
    % Parameters:
    % - img: the image to display (2D array)
    % - spec: the wavelet coefficients (4D array)
    % - Scales: array of scales used in wavelet transform
    % - scale_idx: index of the scale to use
    % - angle_idx: index of the angle to use
    % - clevfactor: factor to adjust contour levels
    % - ProcessFlag: if 1, preprocessing is applied
    
    % Preprocess the image for display
    if ProcessFlag == 1
        img_display = preprocess_img(img);
    else
        img_display = img;
        % Normalize image to [0,1] if not already
        img_display = img_display - min(img_display(:));
        img_display = img_display / max(img_display(:));
    end

    % Display the image
    imagesc(img_display);
    colormap(gray);
    axis image off;  % Set axis aspect ratio and remove axis ticks
    hold on;

    % Extract the wavelet coefficients at the specified scale and angle
    wavelet_real = real(spec(:, :, scale_idx, angle_idx));
    wavelet_abs = abs(spec(:, :, scale_idx, angle_idx));

    % Determine contour levels
    % Use the standard deviation to set contour levels for robustness
    std_real = std(wavelet_real(:));
    std_abs = std(wavelet_abs(:));

    % Set contour levels for real part (crests and troughs)
    num_levels = 5;  % Number of contour levels
    pos_levels = linspace(std_real, max(wavelet_real(:)), num_levels);
    neg_levels = linspace(-std_real, min(wavelet_real(:)), num_levels);

    % Adjust levels with clevfactor
    pos_levels = pos_levels / clevfactor;
    neg_levels = neg_levels / clevfactor;

    % Plot contours of the positive real part (crests)
    contour(wavelet_real, 'LevelList', pos_levels, 'LineColor', 'red', 'LineWidth', 1);

    % Plot contours of the negative real part (troughs)
    contour(wavelet_real, 'LevelList', neg_levels, 'LineColor', 'blue', 'LineWidth', 1);

    % Plot contours of the wavelet power (magnitude squared)
    power_levels = linspace(std_abs^2, max(wavelet_abs(:))^2, num_levels);
    power_levels = power_levels / clevfactor^2;  % Adjust with clevfactor squared

    contour(wavelet_abs.^2, 'LevelList', power_levels, 'LineColor', 'white', 'LineWidth', 1);

    hold off;
end

function data_out = preprocessData(data_in)
    % Clip extremes and apply mild smoothing
    lowerClip = prctile(data_in(:), 1);
    upperClip = prctile(data_in(:), 99);
    data_in(data_in < lowerClip) = lowerClip;
    data_in(data_in > upperClip) = upperClip;

    data_smooth = imgaussfilt(data_in, 1);
    data_out = data_smooth;
end

function img_processed = processDataMethod(data, methodName)
    switch methodName
        case 'raw_normalized'
            img_processed = normalizeData(data);
            img_processed = 1 - img_processed; % invert
            img_processed = img_processed'; % single-channel, transposed as before

        case 'truncated'
            lower_bound = 280; 
            upper_bound = 292.5;
            img_p = data;
            img_p(img_p < lower_bound) = lower_bound;
            img_p(img_p > upper_bound) = upper_bound;
            img_p = (img_p - lower_bound) / (upper_bound - lower_bound);
            img_p = 1 - img_p;
            img_processed = img_p';

        case 'highpass_20'
            img_p = applyHighPass(data, 20, false);
            img_processed = img_p';

        case 'highpass_100'
            img_p = applyHighPass(data, 100, false);
            img_processed = img_p';

        case 'highpass_100_sqrt'
            img_p = applyHighPass(data, 100, true);
            img_processed = img_p';

        case 'highpass_50'
            img_p = applyHighPass(data, 50, false);
            img_processed = img_p';

        case 'highpass_50_sqrt'
            img_p = applyHighPass(data, 50, true);
            img_processed = img_p';

        case 'highpass_20_sqrt'
            img_p = applyHighPass(data, 20, true);
            img_processed = img_p';

        case 'raw_masked'
            maskLimit = 295;
            data(data > maskLimit) = NaN;
            img_p = normalizeDataNaN(data);
            img_p = 1 - img_p;
            img_processed = img_p';

        case 'raw'
            img_p = normalizeDataNaN(data);
            img_p = 1 - img_p;
            img_processed = img_p';
            
        otherwise
            error('Unknown method: %s', methodName);
    end
end

function img_processed = normalizeData(data)
    dmin = min(data(:));
    dmax = max(data(:));
    img_processed = (data - dmin) / (dmax - dmin);
    img_processed(isnan(img_processed)) = 0;
end

function img_processed = normalizeDataNaN(data)
    dmin = nanmin(data(:));
    dmax = nanmax(data(:));
    img_processed = (data - dmin) / (dmax - dmin);
    img_processed(isnan(img_processed)) = 0;
end

function img_processed = applyHighPass(data, filterWidth, doSqrtEnhance)
    lowPass = imgaussfilt(data, filterWidth);
    highPass = data - lowPass;

    if doSqrtEnhance
        highPass = sqrt(abs(highPass)) .* sign(highPass);
    end

    clipMin = -3; 
    clipMax = 3;  
    highPass(highPass < clipMin) = clipMin;
    highPass(highPass > clipMax) = clipMax;

    img_processed = (highPass - clipMin) / (clipMax - clipMin);
    img_processed = 1 - img_processed;
end