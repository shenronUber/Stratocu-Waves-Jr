% --- Code to analyze data and visualize displacements ---

peak_data = peak_data(~any(isnan(peak_data), 2), :);

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

%% Parameters and Setup (First Part)

startDate = datetime(2023, 10, 11);
endDate = datetime(2023, 10, 14);
downloadDir = 'downloaded_data/IR';
methodName = 'highpass_20_sqrt'; % e.g., 'raw_normalized'

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

% Process First Frame (IR Data -> Frame1)
fileName = fileNames{1};  % Only process the first file
filePath = fullfile(downloadDir, fileName);
fprintf('Processing %s with method %s...\n', fileName, methodName);

try
    data = double(ncread(filePath, variableName));
    data_preprocessed = preprocessData(data);
    Frame1 = processDataMethod(data_preprocessed, methodName); % Single-channel 2D array
    fprintf('Frame1 processed successfully using method: %s\n', methodName);
catch ME
    error('Error processing %s: %s\n', fileName, ME.message);
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
square_size_deg = 10; 
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

% Transpose the frame for portrait orientation
frame1 = Frame1;

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

frame1_display = frame1;

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

%%
% --- Code pour le mode basé sur la longueur d’onde ---

% Suppose que filtered_peak_data est déjà peak_data avec valeurs de cohérence filtrées
filtered_peak_data = peak_data;
filtered_peak_data(:, 7) = filtered_coherence_values;

% Garder uniquement les valeurs valides (cohérences non-NaN)
valid_indices = ~isnan(filtered_coherence_values);
filtered_peak_data = filtered_peak_data(valid_indices, :);

% Définir la longueur d’onde désirée (en km), par ex. 50 km
desired_wavelength = 50;

% On considère qu’il y a de légères variations, donc on cherche la longueur d’onde la plus proche
% On calcule la différence absolue par rapport à desired_wavelength
wavelength_diff = abs(filtered_peak_data(:, 8) - desired_wavelength);

% Trouver l’indice du pic qui a la longueur d’onde la plus proche de desired_wavelength
[~, closest_idx] = min(wavelength_diff);
chosen_wavelength = filtered_peak_data(closest_idx, 8);

% Filtrer peak_data pour ne garder que les lignes ayant cette longueur d’onde (exacte ou très proche)
% Si les longueurs d’onde sont discrètes, on peut considérer equality ou une faible tolérance
tolerance = 0.5; % par exemple ±0.5 km
valid_indices = abs(filtered_peak_data(:, 8) - chosen_wavelength) <= tolerance;
desired_wavelength_peaks = filtered_peak_data(valid_indices, :);

% Nombre de frames et squares
num_frames = max(desired_wavelength_peaks(:, 1));
num_squares = max(desired_wavelength_peaks(:, 2));

% Remplir peaks_by_square
peaks_by_square = cell(num_squares, 1);
for s = 1:num_squares
    peaks_by_square{s} = cell(num_frames, 1);
end

for i = 1:size(desired_wavelength_peaks, 1)
    f_idx = desired_wavelength_peaks(i, 1);
    s_idx = desired_wavelength_peaks(i, 2);
    peaks_by_square{s_idx}{f_idx} = desired_wavelength_peaks(i, :);
end

% Calcul des déplacements entre frames
displacements = [];
for s = 1:num_squares
    square_peaks = peaks_by_square{s};
    for frame_idx = 1:(num_frames - 1)
        peak_current = square_peaks{frame_idx};
        peak_next = square_peaks{frame_idx + 1};
        if ~isempty(peak_current) && ~isempty(peak_next)
            angle_current = peak_current(4);
            angle_next = peak_next(4);
            angle_diff = atan2(sin(angle_next - angle_current), cos(angle_next - angle_current));

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
            displacement.speed_current = peak_current(6);
            displacement.speed_next = peak_next(6);

            displacements = [displacements; displacement];
        end
    end
end

% Filtrage des déplacements significatifs
angle_threshold = pi / 8;
wavelength_threshold = 0.1 * mean([displacements.wavelength_current]);
significant_displacements = displacements(abs([displacements.angle_diff]) <= angle_threshold & ...
                                          [displacements.wavelength_diff] <= wavelength_threshold);

% Pour chaque square, on détermine le mode d’angle et la vitesse moyenne associée
frame1_display = frame1;
frame1_display = apply_radial_window(frame1_display, radius_factor, decay_rate);

figure;
imshow(frame1_display, []);
hold on;

all_squares = unique([significant_displacements.square_idx]);
for s = all_squares
    sq_disp = significant_displacements([significant_displacements.square_idx] == s);

    if isempty(sq_disp)
        continue;
    end

    % On veut le mode d'angle_current
    all_angles = [sq_disp.angle_current];

    % Comme les angles sont déjà discrets dans votre dataset, on peut trouver le mode via "unique" et "histc"
    unique_angles = unique(all_angles);
    angle_counts = histc(all_angles, unique_angles);

    [~, max_idx] = max(angle_counts);
    mode_angle = unique_angles(max_idx);

    % Extraire toutes les lignes avec cet angle
    mode_disp = sq_disp([sq_disp.angle_current] == mode_angle);

    % Calculer la vitesse moyenne
    mean_speed = mean([mode_disp.speed_current]);

    % Afficher la flèche sur l'image
    square = squares(s);
    x_center = mean(square.x_range);
    y_center = mean(square.y_range);

    arrow_length = 20;
    dx = arrow_length * cos(mode_angle);
    dy = arrow_length * sin(mode_angle);

    quiver(x_center, y_center, dx, dy, 'LineWidth', 2, 'Color', 'red', 'MaxHeadSize', 2);
    text(x_center + dx, y_center + dy, sprintf('V: %.1f m/s\nλ: %.1f km', mean_speed, chosen_wavelength), ...
        'Color', 'yellow', 'FontSize', 8);
end

title(sprintf('Dominant Direction per Square for λ ≈ %.1f km', chosen_wavelength));
hold off;
%%

filtered_peak_data = peak_data;
filtered_peak_data(:, 7) = filtered_coherence_values;
valid_indices = ~isnan(filtered_coherence_values);
filtered_peak_data = filtered_peak_data(valid_indices, :);

% Define desired angle in DEGREES (e.g., 45°)
desired_angle_deg = 45;

% Convert desired angle to radians
desired_angle = deg2rad(desired_angle_deg);

% Find the dataset angle closest to desired_angle
all_angles = filtered_peak_data(:, 4); % Angles in radians in peak_data
[~, closest_idx] = min(abs(all_angles - desired_angle));
chosen_angle = all_angles(closest_idx); % in radians

% Filter to keep only peaks with angles within a given tolerance
angle_tolerance_deg = 5; % e.g., ±5°
angle_tolerance = deg2rad(angle_tolerance_deg);

valid_indices = abs(filtered_peak_data(:, 4) - chosen_angle) <= angle_tolerance;
desired_angle_peaks = filtered_peak_data(valid_indices, :);

num_frames = max(desired_angle_peaks(:, 1));
num_squares = max(desired_angle_peaks(:, 2));

peaks_by_square = cell(num_squares, 1);
for s = 1:num_squares
    peaks_by_square{s} = cell(num_frames, 1);
end

for i = 1:size(desired_angle_peaks, 1)
    f_idx = desired_angle_peaks(i, 1);
    s_idx = desired_angle_peaks(i, 2);
    peaks_by_square{s_idx}{f_idx} = desired_angle_peaks(i, :);
end

% Calculate displacements
displacements = [];
for s = 1:num_squares
    square_peaks = peaks_by_square{s};
    for frame_idx = 1:(num_frames - 1)
        peak_current = square_peaks{frame_idx};
        peak_next = square_peaks{frame_idx + 1};
        if ~isempty(peak_current) && ~isempty(peak_next)
            angle_current = peak_current(4);
            angle_next = peak_next(4);
            angle_diff = atan2(sin(angle_next - angle_current), cos(angle_next - angle_current));

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
            displacement.speed_current = peak_current(6);
            displacement.speed_next = peak_next(6);

            displacements = [displacements; displacement];
        end
    end
end

% Significant displacement filtering
angle_threshold = deg2rad(5); 
wavelength_threshold = 0.1 * mean([displacements.wavelength_current]);
significant_displacements = displacements(abs([displacements.angle_diff]) <= angle_threshold & ...
                                          [displacements.wavelength_diff] <= wavelength_threshold);

% For each square, find the mode wavelength at this angle and compute mean speed
frame1_display = frame1;
frame1_display = apply_radial_window(frame1_display, radius_factor, decay_rate);

figure;
imshow(frame1_display, []);
hold on;

all_squares = unique([significant_displacements.square_idx]);
for s = all_squares
    sq_disp = significant_displacements([significant_displacements.square_idx] == s);

    if isempty(sq_disp)
        continue;
    end

    % All displacements are around the chosen_angle
    all_wavelengths = [sq_disp.wavelength_current];
    unique_wv = unique(all_wavelengths);
    wv_counts = histc(all_wavelengths, unique_wv);

    [~, max_idx] = max(wv_counts);
    mode_wavelength = unique_wv(max_idx);

    % Extract lines with this mode wavelength
    mode_wv_disp = sq_disp([sq_disp.wavelength_current] == mode_wavelength);

    % Mean speed
    mean_speed = mean([mode_wv_disp.speed_current]);

    % If mean_speed is negative, reverse angle by adding pi (180°) and take absolute value of speed
    display_angle = chosen_angle; % Start with the chosen_angle
    if mean_speed < 0
        mean_speed = abs(mean_speed);
        display_angle = display_angle + pi;
    end

    % Display the arrow: anticlockwise convention
    square = squares(s);
    x_center = mean(square.x_range);
    y_center = mean(square.y_range);

    arrow_length = 20;
    dx = arrow_length * cos(display_angle);
    dy = -arrow_length * sin(display_angle); % Negative sign for anticlockwise, standard mathematical orientation

    quiver(x_center, y_center, dx, dy, 'LineWidth', 2, 'Color', 'red', 'MaxHeadSize', 2);
    text(x_center + dx, y_center + dy, sprintf('V: %.1f m/s\nλ: %.1f km', mean_speed, mode_wavelength), ...
        'Color', 'yellow', 'FontSize', 8);
end

title(sprintf('Dominant Wavelength per Square for Angle ≈ %.1f° (Anticlockwise)', rad2deg(chosen_angle)));
hold off;

%% Vidéo annotée à partir d'une vidéo .mj2 existante et de données peak_data
% Annotates only peaks at the closest wavelength.

% Paramètres
input_video_file = 'background_video.mj2';  % Fichier vidéo d'entrée
output_video_file = 'annotated_video2.mj2';  % Fichier vidéo de sortie annoté
arrow_length = 100;
angle_correction_if_negative_speed = pi; % 180° en radians

% Longueur d'onde désirée (en km)
desired_wavelength = 95; 

% Trouver la longueur d'onde la plus proche dans peak_data
unique_wavelengths = unique(peak_data(:, 8));
[~, closest_idx] = min(abs(unique_wavelengths - desired_wavelength));
closest_wavelength = unique_wavelengths(closest_idx);

disp(['Closest wavelength to ', num2str(desired_wavelength), ' km is ', num2str(closest_wavelength), ' km']);

% Lecture de la vidéo d'entrée
videoReader = VideoReader(input_video_file);

% Création du writer pour la vidéo annotée
videoWriter = VideoWriter(output_video_file, 'Motion JPEG 2000');
videoWriter.FrameRate = videoReader.FrameRate;
videoWriter.LosslessCompression = true;
open(videoWriter);

% Détermination du nombre de frames et de carrés dans peak_data
num_frames = max(peak_data(:,1));
num_squares = max(peak_data(:,2));

frame_idx = 0;
while hasFrame(videoReader)
    frame = readFrame(videoReader);
    frame_idx = frame_idx + 1;

    % Optionnel: s'arrêter un frame avant la fin
    if frame_idx > num_frames - 1
        break; 
    end

    % Création d'une figure temporaire
    h_fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 1920 1080]);
    imshow(frame, []);
    hold on;

    % Annotation de chaque carré
    for s = 1:num_squares
        % Sélection des lignes correspondant à ce frame et ce carré, et à la longueur d'onde la plus proche
        row_mask = (peak_data(:,1) == frame_idx) & (peak_data(:,2) == s) & ...
                   (peak_data(:,8) == closest_wavelength);

        if ~any(row_mask)
            % Aucune donnée pour ce carré à ce frame
            continue;
        end

        % Extraire toutes les lignes correspondantes
        peaks = peak_data(row_mask, :);

        for row_n = 1:size(peaks,1)
            angle = peaks(row_n,4);       % Angle en radians
            speed = peaks(row_n,6);       % Vitesse (peut être négative)
            wavelength = peaks(row_n,8);  % Longueur d'onde (km)
            coherence = peaks(row_n,7);   % Coherence si besoin

            % Si vitesse négative, corriger l'angle et prendre la valeur absolue de la vitesse
            if speed < 0
                speed = abs(speed);
                angle = angle + angle_correction_if_negative_speed;
            end

            % Coordonnées du carré
            square_info = squares(s);
            x_center = mean(square_info.x_range);
            y_center = mean(square_info.y_range);

            % Rééchantillonnage si nécessaire
            x_center_original = x_center * shrinkfactor;
            y_center_original = y_center * shrinkfactor;

            % Calcul de la direction de la flèche
            dx = arrow_length * cos(angle);
            dy = arrow_length * sin(angle);

            % Dessiner la flèche
            quiver(x_center_original, y_center_original, dx, dy, 'LineWidth', 5, 'Color', 'red', 'MaxHeadSize', 2);

            % Ajouter l'annotation textuelle
            text(x_center_original + dx, y_center_original + dy, ...
                sprintf('V: %.1f m/s\n\\lambda: %.1f km', speed, wavelength), ...
                'Color', 'yellow', 'FontSize', 14, 'Interpreter', 'tex');
        end
    end

    title(sprintf('Frame %d | λ ≈ %.1f km', frame_idx, closest_wavelength));
    hold off;

    % Capturer le frame annoté et l'enregistrer
    frame_for_video = getframe(h_fig);
    writeVideo(videoWriter, frame_for_video);

    % Fermer la figure temporaire
    close(h_fig);
end

% Fermer la vidéo annotée
close(videoWriter);

disp('Annotated video created successfully with closest wavelength filtering.');

%% Annotation vidéo basée sur le filtrage par angle

% Paramètres
input_video_file = 'background_video.mj2';      % Vidéo d'entrée
output_video_file = 'annotated_video.mj2';      % Vidéo annotée en sortie
desired_angle_deg = 45;                         % Angle désiré en degrés (sens antihoraire)
arrow_length = 100;
angle_correction_if_negative_speed = pi;        % 180° en radians

% Conversion de l'angle désiré en radians
desired_angle = deg2rad(desired_angle_deg);

% Trouver l'angle le plus proche dans peak_data
unique_angles = unique(peak_data(:, 4));
[~, closest_idx] = min(abs(unique_angles - desired_angle));
chosen_angle = unique_angles(closest_idx);

% Lecture de la vidéo d'entrée
videoReader = VideoReader(input_video_file);

% Création du writer pour la vidéo annotée
videoWriter = VideoWriter(output_video_file, 'Motion JPEG 2000');
videoWriter.FrameRate = videoReader.FrameRate;
videoWriter.LosslessCompression = true;
open(videoWriter);

% Détermination du nombre de frames et de carrés
num_frames = max(peak_data(:,1));
num_squares = max(peak_data(:,2));

frame_idx = 0;
while hasFrame(videoReader)
    frame = readFrame(videoReader);
    frame_idx = frame_idx + 1;

    % Optionnel: s'arrêter un frame avant la fin
    if frame_idx > num_frames - 1
        break; 
    end

    % Création d'une figure temporaire
    h_fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 1920 1080]);
    imshow(frame, []);
    hold on;

    % Annotation de chaque carré
    for s = 1:num_squares
        % Sélection des lignes correspondant à ce frame, ce carré et l'angle choisi
        row_mask = (peak_data(:,1) == frame_idx) & (peak_data(:,2) == s) & ...
                   (peak_data(:,4) == chosen_angle);

        if ~any(row_mask)
            continue;
        end

        % Extraire toutes les lignes correspondantes
        peaks = peak_data(row_mask, :);

        for row_n = 1:size(peaks,1)
            % Correction de l'angle pour passer en sens antihoraire
            angle = -peaks(row_n,4);        % Angle en radians (inversé pour sens antihoraire)
            speed = peaks(row_n,6);        % Vitesse (peut être négative)
            wavelength = peaks(row_n,8);   % Longueur d'onde (km)

            % Si la vitesse est négative, on corrige l'angle
            if speed < 0
                speed = abs(speed);
                angle = angle + angle_correction_if_negative_speed;
            end

            % Coordonnées du carré
            square_info = squares(s);
            x_center = mean(square_info.x_range);
            y_center = mean(square_info.y_range);

            % Rééchantillonnage si nécessaire
            x_center_original = x_center * shrinkfactor;
            y_center_original = y_center * shrinkfactor;

            % Calcul direction flèche
            dx = arrow_length * cos(angle);
            dy = arrow_length * sin(angle);

            % Dessin de la flèche
            quiver(x_center_original, y_center_original, dx, dy, 'LineWidth', 5, 'Color', 'red', 'MaxHeadSize', 2);

            % Position dynamique du label pour éviter les superpositions
            label_offset = 20 * row_n; % Décalage vertical basé sur l'index
            text(x_center_original + dx, y_center_original + dy + label_offset, ...
                sprintf('V: %.1f m/s\n\\lambda: %.1f km', speed, wavelength), ...
                'Color', 'yellow', 'FontSize', 14, 'Interpreter', 'tex');
        end
    end

    % Ajouter un titre pour chaque frame
    title(sprintf('Frame %d | Angle ≈ %.1f°', frame_idx, rad2deg(chosen_angle)));
    hold off;

    % Capturer le frame annoté et l'enregistrer
    frame_for_video = getframe(h_fig);
    writeVideo(videoWriter, frame_for_video);

    % Fermer la figure temporaire
    close(h_fig);
end

% Fermer la vidéo annotée
close(videoWriter);

disp('Annotated video created successfully with single best-matching angle.');


%% Functions

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