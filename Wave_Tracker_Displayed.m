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
        end
    
        threshold = prctile(data(:), 5);
        data(data < threshold) = NaN;
        nan_mask = isnan(data);
        data(isnan(data)) = mean(data,'all','omitnan');
    
        % Processing chain
        img_processed = processDataMethod(data, methodName);
        img_processed(nan_mask')=prctile(img_processed(:), 50);

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
% outputVideoFile = 'all_frames_with_grid.mj2'; 
% v = VideoWriter(outputVideoFile, 'Motion JPEG 2000');
% v.FrameRate = 2;  
% v.LosslessCompression = true; 
% open(v);
% 
% figVid = figure('Position', [100, 100, 800, 600]); 
% 
% for f_idx = 1:num_frames
%     img = processedFrames{f_idx};
%     if shrinkfactor ~= 1
%         img = imresize(img, invshrinkfactor);
%     end
%     if Preprocess_Flag
%         img = preprocess_img(img);
%     end
%     img_windowed = apply_radial_window(img, radius_factor, decay_rate);
%     img_windowed = img_windowed(y_buffer_range, x_buffer_range);
% 
%     imshow(img_windowed, [0, 1], 'Parent', gca);
%     colormap('gray');
%     hold on;
% 
%     for s_idx = 1:length(squares)
%         rectangle('Position', [squares(s_idx).x_range(1) - window_buffer, ...
%                                squares(s_idx).y_range(1) - window_buffer, ...
%                                length(squares(s_idx).x_range), ...
%                                length(squares(s_idx).y_range)], ...
%                   'EdgeColor', 'r', 'LineWidth', 1);
%         text(double(squares(s_idx).x_range(1) - window_buffer), ...
%              double(squares(s_idx).y_range(1) - window_buffer - 5), ...
%              sprintf('%d', squares(s_idx).index), 'Color', 'yellow', 'FontSize', 8);
%     end
% 
%     frame_timestamp = fileTimestamps(f_idx);
%     title(sprintf('Frame %d - %s', f_idx, datestr(frame_timestamp, 'yyyy-mm-dd HH:MM:SS')));
% 
%     hold off;
%     frame_img = getframe(figVid);
%     writeVideo(v, frame_img);
% end
% 
% close(v);
% fprintf('Video created: %s\n', outputVideoFile);

%% ========================================================================
%  Paramètres & Données d'entrée
%  ========================================================================

% Exemple : on va analyser la paire (f, f+1)
f = 12;  % Frame de départ
% Wavelength cible
targetWavelengthKm = 100;  
% Seuil relatif sur la cohérence (type 0.6 => 60% max)
coherenceThreshold = 0.6; 
% Squares d'intérêt pour superposer des contours
squaresOfInterest = [15, 16, 21, 22];

% Fichier .mat où se trouve cwtFullDataMerged
waveletMatFileName = 'allSpecData.mat';

% Ouvre le fichier de CWT
mf = matfile(waveletMatFileName, 'Writable', false);

%% ========================================================================
%  1) Vérifications sur f et squaresOfInterest
%  ========================================================================
num_frames = length(processedFrames);
if f < 1 || (f+1) > num_frames
    error('La paire (f=%d, f+1=%d) est hors bornes (1..%d).', f, f+1, num_frames);
end

%% =========================================================================
%  2) Préparer la frame f pour l'affichage global
%  =========================================================================
frame1 = processedFrames{f};
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

figure;
imshow(frame1_windowed, [0, 1]);
colormap('gray');

hold on;
for idx = 1:length(squares)
    rectangle('Position', [squares(idx).x_range(1), squares(idx).y_range(1), ...
        length(squares(idx).x_range), length(squares(idx).y_range)], ...
        'EdgeColor', 'r', 'LineWidth', 1);
    text(squares(idx).x_range(1), squares(idx).y_range(1) - 5, ...
        sprintf('%d', squares(idx).index), 'Color', 'yellow', 'FontSize', 8);
end
title(sprintf('Frame %d with squares + WaveBlobs', f));

%% =========================================================================
%  3) Charger la CWT (frame f et f+1)
%     On ne fait pas la cohérence "globale", on ira localement square par square
%  =========================================================================
spec1_full_cell = mf.cwtFullDataMerged(f, 1);
spec1_full = spec1_full_cell{1}; % [Ny, Nx, Nscales, Nangles]
spec2_full_cell = mf.cwtFullDataMerged(f+1, 1);
spec2_full = spec2_full_cell{1};

%% =========================================================================
%  4) Boucle: pour chaque square d'intérêt, on calcule la cohérence LOCALE
%     et on utilise find_peaks_and_speeds, puis on retient le pic
%     dont la scale est la plus proche de la targetWavelengthKm
%  =========================================================================
for sqID = squaresOfInterest
    % Récupère la structure du square
    idxSq = find([squares.index] == sqID, 1);
    if isempty(idxSq)
        warning('Square %d introuvable.', sqID);
        continue;
    end
    sq = squares(idxSq);

    % Extraire la portion (y_range, x_range)
    yR = sq.y_range;
    xR = sq.x_range;

    % Sécuriser si besoin (ne pas dépasser la taille de spec1_full)
    yR(yR<1) = [];  yR(yR>size(spec1_full,1)) = [];
    xR(xR<1) = [];  xR(xR>size(spec1_full,2)) = [];

    if isempty(yR) || isempty(xR)
        warning('Square %d => range out of bounds, skip.', sqID);
        continue;
    end

    % sp1, sp2: local wavelet portion
    sp1 = spec1_full(yR, xR, :, :); % [height, width, Nscales, Nangles]
    sp2 = spec2_full(yR, xR, :, :);

    % xwt => cohérence
    xwt = sp1 .* conj(sp2);
    coherence = abs(xwt);
    phase = angle(xwt);

    % (Optionnel) Normalisation "plot_waverose"-like => on peut la faire
    % scale par scale. Mais si vous voulez EXACTEMENT la logique "plot_waverose",
    % vous pouvez faire:
    for sc = 1:length(Scales)
        tmp = coherence(:,:,sc,:);
        scmean = mean(tmp(:), 'omitnan');
        if scmean ~= 0
            coherence(:,:,sc,:) = coherence(:,:,sc,:) / scmean;
        end
    end

    % On appelle find_peaks_and_speeds => (Scales, Angles, pixel_size_km, time_interval)
    peak_list = find_peaks_and_speeds(coherence, phase, ...
        Scales, Angles, pixel_size_km, time_interval);
    % peak_list est typiquement un tableau (Nscales x 6), 
    %  colonnes = [ scale_px, angle_rad, mean_phase, speed_m_s, max_coh_value, wavelength_km ]

    % Chercher la ligne dont la "wavelength_km" est la plus proche de targetWavelengthKm
    % Le 6ème champ de peak_list(:,:) = wavelength_km
    diffs = abs(peak_list(:,6) - targetWavelengthKm);
    [~, bestIdx] = min(diffs);
    chosenScale = peak_list(bestIdx, 1);  % scale en px
    chosenAngle = peak_list(bestIdx, 2);
    chosenWave  = peak_list(bestIdx, 6);

    fprintf('\nSquare #%d => bestPeakIdx=%d, scale=%.2f px (%.1f km), angle=%.2f rad\n', ...
        sqID, bestIdx, chosenScale, chosenWave, chosenAngle);

    % On retrouve l'indice scale_idx & angle_idx
    %  => scale_idx = l'indice i tel que Scales(i) ~ chosenScale
    [~, scale_idx] = min(abs(Scales - chosenScale));
    chosenAngle = mod(chosenAngle, 2*pi);
    [~, angle_idx] = min(abs(Angles - chosenAngle));

    % Extraire la cohSlice correspondante
    cohSlice = coherence(:,:,scale_idx, angle_idx); % 2D
    cmax = max(cohSlice(:));
    if cmax <= 0
        fprintf('  => cmax=%.3g => skip.\n', cmax);
        continue;
    end

    % On fabrique un masque binaire => coherenceThreshold * cmax
    thrVal = coherenceThreshold * cmax;
    mask = (cohSlice >= thrVal);

    % => On veut tracer le contour de ce masque en magenta
    %    Mais sur la figure "globale"
    % => Le square est aux coords [xR, yR], or le mask = [size(yR), size(xR)]
    % => On utilise contour() avec X=meshgrid(xR, yR), Z=cohSlice
    [Xgrid, Ygrid] = meshgrid(xR, yR);

    % On ne fait pas contour(cohSlice, [thrVal thrVal]) => on fait:
    contour(Xgrid, Ygrid, cohSlice, [thrVal thrVal], 'm', 'LineWidth', 2);

    % Optionnel: marquer centroids
    CC = bwconncomp(mask);
    S = regionprops(CC, 'Centroid');
    for b = 1:length(S)
        cLocal = S(b).Centroid;  % [x, y] dans "cohSlice"
        % Convertir en coords globales
        real_x = xR(1) + cLocal(1) - 1;
        real_y = yR(1) + cLocal(2) - 1;
        plot(real_x, real_y, 'mo', 'MarkerSize', 3, 'LineWidth', 1.5);
    end

    fprintf('  => cmax=%.4f, threshold=%.4f => drawn magenta contour.\n', cmax, thrVal);
end

hold off;

%%

% Nom du fichier vidéo de sortie
outputVideoFile = 'wavelet_analysis_video.avi';

% Initialiser l'objet VideoWriter
videoObj = VideoWriter(outputVideoFile, 'Motion JPEG AVI');
videoObj.FrameRate = 5;  % Ajustez la vitesse (frames/sec)
open(videoObj);

% Définir la plage de frames à analyser
frameRange = 1:(num_frames - 1);  % Paires (f, f+1)

% Boucle sur les frames
for f = frameRange
    % Créer une nouvelle figure
    figure('Visible', 'off', 'Color', 'w');
    
    % Préparer la frame (f) pour l'affichage
    frame1 = processedFrames{f};
    if shrinkfactor ~= 1
        frame1 = imresize(frame1, invshrinkfactor);
    end
    if Preprocess_Flag
        frame1 = preprocess_img(frame1);
    end
    frame1_windowed = apply_radial_window(frame1, radius_factor, decay_rate);
    
    % Afficher la frame (avec les squares en rouge)
    imshow(frame1_windowed, [0, 1]);
    colormap('gray');
    hold on;
    for idx = 1:length(squares)
        rectangle('Position', [squares(idx).x_range(1), squares(idx).y_range(1), ...
            length(squares(idx).x_range), length(squares(idx).y_range)], ...
            'EdgeColor', 'r', 'LineWidth', 1);
        text(squares(idx).x_range(1), squares(idx).y_range(1) - 5, ...
            sprintf('%d', squares(idx).index), 'Color', 'yellow', 'FontSize', 8);
    end
    title(sprintf('Frame %d with squares + WaveBlobs', f));
    
    % Charger les CWT (frames f et f+1)
    spec1_full_cell = mf.cwtFullDataMerged(f, 1);
    spec1_full = spec1_full_cell{1}; % [Ny, Nx, Nscales, Nangles]
    spec2_full_cell = mf.cwtFullDataMerged(f+1, 1);
    spec2_full = spec2_full_cell{1};
    
    % Boucle sur les squares d'intérêt
    for sqID = squaresOfInterest
        % Récupérer les coordonnées du square
        idxSq = find([squares.index] == sqID, 1);
        if isempty(idxSq)
            warning('Square %d introuvable.', sqID);
            continue;
        end
        sq = squares(idxSq);

        % Récupérer la portion (y_range, x_range)
        yR = sq.y_range;
        xR = sq.x_range;
        yR(yR<1) = [];  yR(yR>size(spec1_full,1)) = [];
        xR(xR<1) = [];  xR(xR>size(spec1_full,2)) = [];
        if isempty(yR) || isempty(xR), continue; end

        % Extraire les portions locales sp1 et sp2
        sp1 = spec1_full(yR, xR, :, :);
        sp2 = spec2_full(yR, xR, :, :);

        % Calculer cohérence et phase
        xwt = sp1 .* conj(sp2);
        coherence = abs(xwt);
        phase = angle(xwt);

        % Normaliser par échelle
        for sc = 1:length(Scales)
            tmp = coherence(:,:,sc,:);
            scmean = mean(tmp(:), 'omitnan');
            if scmean ~= 0
                coherence(:,:,sc,:) = coherence(:,:,sc,:) / scmean;
            end
        end

        % Appeler find_peaks_and_speeds
        peak_list = find_peaks_and_speeds(coherence, phase, ...
            Scales, Angles, pixel_size_km, time_interval);

        % Chercher le pic le plus proche de targetWavelengthKm
        diffs = abs(peak_list(:,6) - targetWavelengthKm);
        [~, bestIdx] = min(diffs);
        chosenScale = peak_list(bestIdx, 1);  % scale en px
        chosenAngle = peak_list(bestIdx, 2);

        % Identifier les indices de scale et angle
        [~, scale_idx] = min(abs(Scales - chosenScale));
        chosenAngle = mod(chosenAngle, 2*pi);
        [~, angle_idx] = min(abs(Angles - chosenAngle));

        % Extraire cohSlice
        cohSlice = coherence(:,:,scale_idx, angle_idx);
        cmax = max(cohSlice(:));
        if cmax <= 0, continue; end

        % Masque pour le seuil
        thrVal = coherenceThreshold * cmax;
        mask = (cohSlice >= thrVal);

        % Tracer le contour sur la figure
        [Xgrid, Ygrid] = meshgrid(xR, yR);
        contour(Xgrid, Ygrid, cohSlice, [thrVal thrVal], 'm', 'LineWidth', 2);
    end
    hold off;

    % Capturer la frame pour la vidéo
    frame = getframe(gcf);
    writeVideo(videoObj, frame);

    % Fermer la figure courante
    close(gcf);
end

% Fermer l'objet vidéo
close(videoObj);
disp(['Vidéo sauvegardée dans ', outputVideoFile]);



%% Functions

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

function img_processed = preprocess_img(img)
    lower_bound = prctile(img(:), 10);
    upper_bound = prctile(img(:), 99);
    img_processed = img;
    img_processed(img < lower_bound) = lower_bound;
    img_processed(img > upper_bound) = upper_bound;
    img_processed = (img_processed - lower_bound) / (upper_bound - lower_bound);
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
