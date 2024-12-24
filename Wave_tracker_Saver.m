%% Parameters and Setup (First Part)

startDate = datetime(2023, 10, 12, 1, 0, 0);
endDate   = datetime(2023, 10, 12, 9, 0, 0);
downloadDir = 'downloaded_data/IR';
methodName = 'highpass_50_sqrt'; % e.g., 'raw_normalized' or 'highpass_20_sqrt'

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
        
        % Replace NaNs with zero and print how many
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
total_time = tic;

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

% Zone d'intérêt pour une résolution accrue
focus_min = 40;   % Échelle minimum d'intérêt
focus_max = 125;  % Échelle maximum d'intérêt
n_focus_scales = 10; % Nombre d'échelles dans la zone d'intérêt
n_outer_scales = 5+1;  % Nombre d'échelles dans les zones extérieures

% Génération des échelles
outer_scales_low  = logspace(log10(min_scale_km),  log10(focus_min),  n_outer_scales);
focus_scales      = linspace(focus_min, focus_max, n_focus_scales);
outer_scales_high = logspace(log10(focus_max), log10(max_scale_km),  n_outer_scales);

% Suppression des doublons aux jonctions
outer_scales_low = outer_scales_low(1:end-1);  % Exclure la dernière valeur
outer_scales_high = outer_scales_high(2:end);  % Exclure la première valeur

% Nombre total d'échelles (approx.)
Scales_km = [outer_scales_low, focus_scales, outer_scales_high];
NSCALES = numel(Scales_km);
Scales = Scales_km / pixel_size_km; % Conversion en pixels

window_buffer = 10; 
Preprocess_Flag = 1; 
time_interval = 1800; 
Save_Metadata_Flag = 0; % si on veut stocker des données additionnelles dans 'square_metadata'

% Radial window parameters
radius_factor = 1.1;  
decay_rate = 0.05;   

%% Prepare First Frame for Dimension Calculation
frame1 = processedFrames{1};
if shrinkfactor ~= 1
    frame1 = imresize(frame1, invshrinkfactor);
end
if Preprocess_Flag
    frame1 = preprocess_img(frame1);
end
frame1_windowed = apply_radial_window(frame1, radius_factor, decay_rate);

[frame_height_final, frame_width_final] = size(frame1_windowed);
x_buffer_range = (window_buffer + 1):(frame_width_final - window_buffer);
y_buffer_range = (window_buffer + 1):(frame_height_final - window_buffer);

adjusted_frame_width = length(x_buffer_range);
adjusted_frame_height = length(y_buffer_range);

square_size_px = round(square_size_km / pixel_size_km);
num_squares_x = ceil(adjusted_frame_width / square_size_px);
num_squares_y = ceil(adjusted_frame_height / square_size_px);

% Generate Squares
squares = [];
idx = 1;
for i = 1:num_squares_y
    for j = 1:num_squares_x
        x_start = floor((j - 1) * adjusted_frame_width / num_squares_x) + 1;
        y_start = floor((i - 1) * adjusted_frame_height / num_squares_y) + 1;
        x_end   = floor(j * adjusted_frame_width / num_squares_x);
        y_end   = floor(i * adjusted_frame_height / num_squares_y);

        if x_end > x_start && y_end > y_start
            squares(idx).x_range = x_buffer_range(x_start:x_end);
            squares(idx).y_range = y_buffer_range(y_start:y_end);
            squares(idx).index   = idx;
            idx = idx + 1;
        end
    end
end

%% Display Frame 1 with Squares
figure;
if Preprocess_Flag
    frame1 = preprocess_img(frame1); % on refait un preprocess si besoin
end
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

%% Create a video of all frames with the grid overlay
outputVideoFile = 'all_frames_with_grid.mj2'; 
v = VideoWriter(outputVideoFile, 'Motion JPEG 2000');
v.FrameRate = 2;  
v.LosslessCompression = true; 
open(v);

figVid = figure('Position', [100, 100, 800, 600]); 

for f_idx = 1:num_frames
    img = processedFrames{f_idx};
    if shrinkfactor ~= 1
        img = imresize(img, invshrinkfactor);
    end
    if Preprocess_Flag
        img = preprocess_img(img);
    end
    img_windowed = apply_radial_window(img, radius_factor, decay_rate);
    img_windowed = img_windowed(y_buffer_range, x_buffer_range);

    imshow(img_windowed, [0, 1], 'Parent', gca);
    colormap('gray');
    hold on;

    for s_idx = 1:length(squares)
        rectangle('Position', [squares(s_idx).x_range(1) - window_buffer, ...
                               squares(s_idx).y_range(1) - window_buffer, ...
                               length(squares(s_idx).x_range), ...
                               length(squares(s_idx).y_range)], ...
                  'EdgeColor', 'r', 'LineWidth', 1);
        text(double(squares(s_idx).x_range(1) - window_buffer), ...
             double(squares(s_idx).y_range(1) - window_buffer - 5), ...
             sprintf('%d', squares(s_idx).index), 'Color', 'yellow', 'FontSize', 8);
    end

    frame_timestamp = fileTimestamps(f_idx);
    title(sprintf('Frame %d - %s', f_idx, datestr(frame_timestamp, 'yyyy-mm-dd HH:MM:SS')));

    hold off;
    frame_img = getframe(figVid);
    writeVideo(v, frame_img);
end

close(v);
fprintf('Video created: %s\n', outputVideoFile);

%% Étape 1 : Calcul et sauvegarde de la CWT pour chaque frame séparément

fprintf('\n=== Computing wavelet transform for each frame and saving to disk ===\n');
for f_idx = 1:num_frames
    img = processedFrames{f_idx};
    
    if shrinkfactor ~= 1
        img = imresize(img, invshrinkfactor);
    end
    if Preprocess_Flag
        img = preprocess_img(img);
    end
    img = apply_radial_window(img, radius_factor, decay_rate);
    img = img(y_buffer_range, x_buffer_range);
    
    % Calcul de la transformée en ondelettes 2D
    cwt_data = cwtft2(img, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
    spec_full = squeeze(cwt_data.cfs);  

    % Sauvegarde immédiate dans un fichier .mat distinct
    filenameSpec = sprintf('specData_frame_%d.mat', f_idx);
    save(filenameSpec, 'spec_full', '-v7.3');
    
    fprintf('Frame %d/%d wavelet transform computed & saved (%s).\n', ...
        f_idx, num_frames, filenameSpec);
end
fprintf('=== Wavelet transforms saved for all frames ===\n\n');

%% Étape 2 (Optionnel) : Fusionner tous les spec_full en un seul fichier
%   puis supprimer les fichiers intermédiaires

MergeAndDelete = true;

if MergeAndDelete
    fprintf('Merging all specData_frame_*.mat into one file: allSpecData.mat\n');
    
    % 1) Créez une variable "cwtFullDataMerged" vide mais dimensionnée
    cwtFullDataMerged = cell(num_frames, 1);  %#ok<NASGU> 
    %   --> Sauvegardez-la dans un fichier v7.3
    save('allSpecData.mat', 'cwtFullDataMerged', '-v7.3');
    
    % 2) Ouvrez ce fichier en mode "Writable"
    mf = matfile('allSpecData.mat', 'Writable', true);
    
    % 3) Parcourez tous les specData_frame_i.mat
    for f_idx = 1:num_frames
        filenameSpec = sprintf('specData_frame_%d.mat', f_idx);
        tmp = load(filenameSpec, 'spec_full');
        
        % Ecriture DIRECTEMENT dans mf.cwtFullDataMerged{f_idx}
        mf.cwtFullDataMerged(f_idx, 1) = {tmp.spec_full};
    end
    
    fprintf('allSpecData.mat created with all frames.\n');
    
    % 4) Suppression des fichiers intermédiaires
    for f_idx = 1:num_frames
        filenameSpec = sprintf('specData_frame_%d.mat', f_idx);
        if exist(filenameSpec, 'file')
            delete(filenameSpec);
        end
    end
    fprintf('Intermediate specData files deleted.\n');
end


%% Étape 3 : Boucle de traitement en paires de frames (Cross-Wavelet)
%  en rechargeant depuis le fichier 'allSpecData.mat' (pas de surcharge mémoire)

% Dimensions du tableau peak_data
total_entries = (num_frames - 1) * length(squares) * NSCALES;
peak_data = NaN(total_entries, 8);
entry_counter = 0;

% Optionnel : structure pour stocker metadata, si besoin
square_metadata = struct('frame_idx', {}, 'square_idx', {}, 'spec1', {});

fprintf('\n=== Processing frames in pairs (Cross-Wavelet) ===\n');

% Ouvrir le fichier unique en lecture
mf = matfile('allSpecData.mat'); 

for frame_idx = 1:(num_frames - 1)
    fprintf('Processing frame pair %d/%d...\n', frame_idx, num_frames - 1);
    frame_time = tic;
    
    % Charger la CWT du frame_idx et frame_idx+1 PARTIELLEMENT
    spec1_full = mf.cwtFullDataMerged{frame_idx};      % 4D array: (y, x, scale, angle)
    spec2_full = mf.cwtFullDataMerged{frame_idx + 1};  % idem
    
    % Charger l'image fenêtrée (frame_idx) pour les tests brightness
    img1 = processedFrames{frame_idx};
    if shrinkfactor ~= 1
        img1 = imresize(img1, invshrinkfactor);
    end
    if Preprocess_Flag
        img1 = preprocess_img(img1);
    end
    img1 = apply_radial_window(img1, radius_factor, decay_rate);
    img1 = img1(y_buffer_range, x_buffer_range);

    for s_idx = 1:length(squares)
        square = squares(s_idx);
        x_range = square.x_range - window_buffer; 
        y_range = square.y_range - window_buffer;
        
        spec1 = spec1_full(y_range, x_range, :, :);
        spec2 = spec2_full(y_range, x_range, :, :);

        if Save_Metadata_Flag == 1
            metadata_entry.frame_idx = frame_idx;
            metadata_entry.square_idx = s_idx;
            metadata_entry.spec1 = spec1; 
            square_metadata = [square_metadata; metadata_entry]; %#ok<AGROW>
        end

        % Test sur la luminosité
        sub_img1 = img1(y_range, x_range);
        mean_brightness = mean(sub_img1(:));
        std_brightness  = std(sub_img1(:));
        
        if mean_brightness < brightness_threshold || std_brightness > std_threshold
            continue; 
        end
        
        % Cross-wavelet
        xwt = spec1 .* conj(spec2);
        coherence = abs(xwt);
        phase_difference = angle(xwt);

        % Calcul de pics et vitesses
        peak_list = find_peaks_and_speeds(coherence, phase_difference, ...
                                          Scales, Angles, pixel_size_km, time_interval);
        n_scales = size(peak_list, 1);

        data_to_store = [repmat(frame_idx, n_scales, 1), ...
                         repmat(s_idx,     n_scales, 1), ...
                         peak_list];
        
        idx_start = entry_counter + 1;
        idx_end   = entry_counter + n_scales;
        
        % Agrandir peak_data si besoin
        if idx_end > total_entries
            extra_entries = total_entries; 
            peak_data = [peak_data; NaN(extra_entries, 8)]; %#ok<AGROW>
            total_entries = total_entries + extra_entries;
        end
        
        peak_data(idx_start:idx_end, :) = data_to_store;
        entry_counter = idx_end;
    end
    
    fprintf('Time to process frame pair %d: %.2f seconds.\n', frame_idx, toc(frame_time));
end

peak_data = peak_data(1:entry_counter, :);

fprintf('Total execution time: %.2f seconds.\n', toc(total_time));

%% (Optionally) Analyze the data
plot_waverose(12,15,false,4,"allSpecData.mat", processedFrames, squares, ...
    Scales, Angles, pixel_size_km, time_interval, ...
    shrinkfactor, invshrinkfactor, Preprocess_Flag, ...
    radius_factor, decay_rate, ...
    x_buffer_range, y_buffer_range, ...
    brightness_threshold, std_threshold);

%% Supporting Functions
function data = readVariableFromFile(filePath)
    info = ncinfo(filePath);
    variables = {info.Variables.Name};
    if any(strcmp(variables, 'CMI'))
        variableName = 'CMI';
    elseif any(strcmp(variables, 'Rad'))
        variableName = 'Rad';
    else
        warning('No suitable variable found in file %s. Available variables: %s.', ...
                filePath, strjoin(variables, ', '));
        error('No suitable variable to read.');
    end
    data = ncread(filePath, variableName);
end

function img_processed = preprocess_img(img)
    lower_bound = prctile(img(:), 10);
    upper_bound = prctile(img(:), 99);
    img_processed = img;
    img_processed(img < lower_bound) = lower_bound;
    img_processed(img > upper_bound) = upper_bound;
    img_processed = (img_processed - lower_bound) / (upper_bound - lower_bound);
end

function peak_list = find_peaks_and_speeds(coherence, phase_difference, Scales, Angles, pixel_size_km, time_interval)
    num_scales = length(Scales);
    peak_list = zeros(num_scales, 6); 
    for scale_idx = 1:num_scales
        coherence_scale = squeeze(mean(mean(coherence(:, :, scale_idx, :), 1, 'omitnan'), 2, 'omitnan'));
        [max_coherence_value, angle_idx] = max(coherence_scale);
        angle = Angles(angle_idx);

        phase_slice = phase_difference(:, :, scale_idx, angle_idx);
        coherence_slice = coherence(:, :, scale_idx, angle_idx);

        coherence_threshold = 0.5 * max(coherence_slice(:));
        coherence_mask = coherence_slice >= coherence_threshold;
        if any(coherence_mask(:))
            mean_phase_diff = mean(phase_slice(coherence_mask), 'omitnan');
        else
            mean_phase_diff = NaN; 
        end

        wavelength_km = Scales(scale_idx)*pi/sqrt(2)*pixel_size_km;
        distance_shift_km = (mean_phase_diff * wavelength_km) / (2 * pi);
        speed_m_per_s = (distance_shift_km / time_interval) * 1000;

        peak_list(scale_idx, :) = [Scales(scale_idx), angle, mean_phase_diff, speed_m_per_s, max_coherence_value, wavelength_km];
    end
end

function img_windowed = apply_radial_window(img, radius_factor, decay_rate)
    [rows, cols] = size(img);
    center_x = cols / 2;
    center_y = rows / 2;
    [x, y] = meshgrid(1:cols, 1:rows);
    distances = sqrt((x - center_x).^2 + (y - center_y).^2);
    max_distance = radius_factor * min(center_x, center_y);
    window = 1 ./ (1 + exp(decay_rate * (distances - max_distance)));
    img_windowed = img .* window;
end

function plot_waverose(frame_id, square_id, display_sine_wave_plots, figs_per_page, ...
    waveletMatFileName, processedFrames, squares, ...
    Scales, Angles, pixel_size_km, time_interval, ...
    shrinkfactor, invshrinkfactor, Preprocess_Flag, ...
    radius_factor, decay_rate, ...
    x_buffer_range, y_buffer_range, ...
    brightness_threshold, std_threshold)
%PLOT_WAVEROSE Generate an advanced rose plot for a specific frame & square,
% en chargeant les CWT depuis un fichier .mat (allSpecData.mat par ex).
%
% SYNTAX:
%   plot_waverose(frame_id, square_id, display_sine_wave_plots, figs_per_page, ...
%       waveletMatFileName, processedFrames, squares, ...
%       Scales, Angles, pixel_size_km, time_interval, ...
%       shrinkfactor, invshrinkfactor, Preprocess_Flag, ...
%       radius_factor, decay_rate, ...
%       x_buffer_range, y_buffer_range, ...
%       brightness_threshold, std_threshold)
%
% INPUTS:
%   - frame_id, square_id, display_sine_wave_plots, figs_per_page: comme avant
%   - waveletMatFileName : nom du .mat avec 'cwtFullDataMerged'
%   - processedFrames, squares, ...
%   - Scales, Angles, pixel_size_km, time_interval, ...
%   - shrinkfactor, invshrinkfactor, Preprocess_Flag, ...
%   - radius_factor, decay_rate, ...
%   - x_buffer_range, y_buffer_range, ...
%   - brightness_threshold, std_threshold
%
% NOTE:
%   Cette version ne calcule pas la CWT localement. Elle la charge depuis
%   waveletMatFileName, en supposant que cwtFullDataMerged{frame_id} =
%   (Ny, Nx, Nscales, Nangles).
%

%% 0) Gérer les arguments optionnels
if nargin < 3 || isempty(display_sine_wave_plots)
    display_sine_wave_plots = true;
end
if nargin < 4 || isempty(figs_per_page)
    figs_per_page = 20;
end

%% 1) Vérifier la validité de frame_id
num_frames = length(processedFrames);
if frame_id < 1 || frame_id >= num_frames
    error('Invalid frame_id. Must be between 1 and %d.', num_frames - 1);
end

%% 2) Charger les CWT depuis le matfile
mf = matfile(waveletMatFileName, 'Writable', false);
% On suppose qu'il y a cwtFullDataMerged dans ce fichier
spec1_full = mf.cwtFullDataMerged{frame_id};
spec2_full = mf.cwtFullDataMerged{frame_id + 1};
% spec1_full, spec2_full : (Ny, Nx, Nscales, Nangles)

%% 3) Récupérer et prétraiter les frames 2D (pour test luminosité)
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

% Buffer
frame1_windowed = frame1_windowed(y_buffer_range, x_buffer_range);
frame2_windowed = frame2_windowed(y_buffer_range, x_buffer_range);

%% 4) Extraire la zone square_id
idx = find([squares.index] == square_id, 1);
if isempty(idx)
    error('Invalid square_id. Square not found.');
end
square = squares(idx);

x_range = square.x_range;
y_range = square.y_range;

% Attention: vérifier que x_range, y_range sont cohérents avec spec1_full
% (par ex. s'ils ont déjà soustrait window_buffer, etc.)

% Filtrer si besoin
x_range(x_range < 1) = [];
y_range(y_range < 1) = [];
x_range(x_range > size(spec1_full, 2)) = [];
y_range(y_range > size(spec1_full, 1)) = [];

% Extraire la portion
spec1 = spec1_full(y_range, x_range, :, :);
spec2 = spec2_full(y_range, x_range, :, :);

% Test luminosité
img1 = frame1_windowed(y_range, x_range);
mean_brightness = mean(img1(:));
std_brightness  = std(img1(:));

if mean_brightness < brightness_threshold || std_brightness > std_threshold
    error('The selected square does not meet brightness or std thresholds.');
end

%% 5) Cross-Wavelet (coherence, phase)
xwt = spec1 .* conj(spec2);
coherence = abs(xwt);
phase_difference = angle(xwt);

% Normalisations éventuelles
for s_idx = 1:length(Scales)
    this_abs_xwt = coherence(:,:,s_idx,:);
    scale_mean = mean(this_abs_xwt(:), 'omitnan'); 
    if scale_mean ~= 0
        coherence(:,:,s_idx,:) = coherence(:,:,s_idx,:) / scale_mean;
    end
end

% Power spectra
power1 = abs(spec1).^2;
for s_idx = 1:length(Scales)
    this_abs_xwt = power1(:,:,s_idx,:);
    scale_mean = sqrt(mean(this_abs_xwt(:).^2, 'omitnan')); 
    if scale_mean ~= 0
        power1(:,:,s_idx,:) = power1(:,:,s_idx,:) / scale_mean;
    end
end

%% 6) Trouver les pics
peak_list = find_peaks_and_speeds(coherence, phase_difference, ...
    Scales, Angles, pixel_size_km, time_interval);

%% 7) Plot de la rose (idem code original)
% Idem, on reprend tout le code de génération de la rose polaire,
% overlays, etc.

% -----------------------------------------------------------------------
% ... (on recopie la partie "Plot Advanced Rose Plot" etc.)
% -----------------------------------------------------------------------

figure;  % etc.

% [NB: Dans l’exemple ici, je mets juste un commentaire. 
%  Tu peux copier/coller l’intégralité du code de tracé que tu avais,
%  en conservant la logique de "anglespec_power" / "anglespec_coherence" etc.,
%  puisque rien ne change, hormis le fait que la CWT est déjà chargée, 
%  au lieu d’être recalculée localement.]

% -----------------------------------------------------------------------
% Idem pour "Coherence Maps", "Phase Difference Maps", etc.
% ...

%% 8) Calcul & Affichage des Speed (sine waves)
num_peaks = size(peak_list, 1);
speeds = zeros(num_peaks, 1);
for i = 1:num_peaks
    scale_px = peak_list(i, 1);
    mean_phase_diff = peak_list(i, 3);

    % Wavelength en km
    wavelength_km = scale_px * pi/sqrt(2) * pixel_size_km;
    distance_shift_km = mean_phase_diff * wavelength_km / (2*pi);
    speed_km_per_s = distance_shift_km / time_interval;
    speeds(i) = speed_km_per_s;
    
    peak_list(i,4) = speed_km_per_s * 1000; % => m/s

    if display_sine_wave_plots
        % idem code existant
        x_values = linspace(0, wavelength_km, 100);
        sine_wave_original = sin(2*pi*x_values / wavelength_km);
        sine_wave_with_phase = sin(2*pi*x_values / wavelength_km + mean_phase_diff);
        
        figure;
        plot(x_values, sine_wave_original, 'b-', 'LineWidth', 2);
        hold on;
        plot(x_values, sine_wave_with_phase, 'r--', 'LineWidth', 2);
        % ...
        title(sprintf('Peak %d: Speed=%.3f m/s', i, speed_km_per_s*1000));
        % ...
    end
end

end % end of function

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
            img_processed = 1 - img_processed;
            img_processed = img_processed';
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
