function MultiInstrument_CrossWavelet()
    %% 0) Data Formatting

    startDate   = datetime(2023,09,28,14,0,0);
    endDate     = datetime(2023,09,28,20,0,0);
    
    sourceRoot   = 'C:\Users\admin\Box\GOES2go_satellite_downloads';
    outputRoot   = 'C:\Users\admin\Box\GWaves_2023_10_11-14_SEPAC';
    
    % For IR
    renameAndOrganizeFiles('IR', startDate, endDate, sourceRoot, outputRoot);
    
    % For VIS
    renameAndOrganizeFiles('VIS', startDate, endDate, sourceRoot, outputRoot);

    % For WV ...


    %% 1) User Parameters
    %----------------------------------------------------------------------
    startDate   = datetime(2023,09,28,14,0,0);
    endDate     = datetime(2023,09,28,20,0,0);
    
    % List of instruments to process (extendable)
    instrumentList = {'IR','VIS'};
    % e.g., in the future: {'IR','VIS','WV'}

    % Root base directory that contains:
    %    <Instrument>\Data
    % into which files named <Instrument>_YYYY_MM_DD_HH_MM.nc exist
    rootBaseDir = 'C:\Users\admin\Box\GWaves_2023_10_11-14_SEPAC';

    % Wavelet parameters
    NANGLES = 12;
    Angles  = linspace(0, pi, NANGLES);  % same as 0 : pi/(NANGLES-1) : pi

    Scales_km = [10,50,80,110,140,170,200,230,300,405]; % example
    NSCALES   = numel(Scales_km);

    % Spatial scaling
    degrees_per_pixel = 0.04;
    km_per_degree     = 111.32;
    original_px_km    = degrees_per_pixel * km_per_degree;  % ~4.45 km
    shrinkfactor      = 2;
    invshrinkfactor   = 1 / shrinkfactor;
    pixel_size_km     = original_px_km * shrinkfactor;      
    % Convert scales to pixel units (after shrink)
    Scales = Scales_km / pixel_size_km; 

    % Radial window parameters (optional)
    doRadialWindow = true;
    radius_factor  = 1.1; 
    decay_rate     = 0.05;

    % Buffer for dividing into squares
    window_buffer = 10;  % number of pixels to ignore on each border
    square_size_deg       = 5;   % dimension of squares in degrees
    effective_deg_per_px  = degrees_per_pixel; 

    %% 2) Gather files for each instrument
    %----------------------------------------------------------------------
    allInstrumentData = struct(); 
    % Each field will be e.g. allInstrumentData.IR with .Files, .Times, .VarName

    for iInst = 1:numel(instrumentList)
        instr   = instrumentList{iInst};
        dataDir = fullfile(rootBaseDir, upper(instr), 'Data');
        if ~exist(dataDir, 'dir')
            error('Data directory for %s not found: %s', instr, dataDir);
        end

        [fNames, fTimes, varName] = getDateRangeFiles(dataDir, startDate, endDate);
        allInstrumentData.(instr).Files   = fNames;
        allInstrumentData.(instr).Times   = fTimes;
        allInstrumentData.(instr).VarName = varName;
    end

    % Intersect the timestamps among all instruments
    % If you have more than two instruments, you can iteratively intersect them all.
    commonTimes = allInstrumentData.(instrumentList{1}).Times;
    for iInst = 2:numel(instrumentList)
        commonTimes = intersect(commonTimes, allInstrumentData.(instrumentList{iInst}).Times);
    end

    % If no common times, exit
    if isempty(commonTimes)
        fprintf('No common timestamps found among all instruments. Exiting.\n');
        return;
    end

    % Create indexing arrays for each instrument to map commonTimes back to file indices
    indexArrays = struct();
    for iInst = 1:numel(instrumentList)
        instr   = instrumentList{iInst};
        iTimes  = allInstrumentData.(instr).Times;
        [~, idxIn] = intersect(iTimes, commonTimes);
        indexArrays.(instr) = idxIn;  % holds the index in fNames/fTimes
    end

    %% 3) Main Loop over each common time
    %----------------------------------------------------------------------
    numCommon = numel(commonTimes);
    fprintf('Found %d common frames among instruments [%s].\n', ...
        numCommon, strjoin(instrumentList, ','));

    for f_idx = 1:numCommon
        thisTime = commonTimes(f_idx);
        frameDateStr = datestr(thisTime, 'yyyy_mm_dd_HHMMSS');
        fprintf('\n[%d/%d] Processing time = %s\n', f_idx, numCommon, frameDateStr);

        % Step 1: Load data for each instrument and do single wavelet if needed
        waveletSpecs = struct();  % will hold spec_full for each instrument
        squaresInfo  = struct();  % we only need to compute squares once, but let's store them

        for iInst = 1:numel(instrumentList)
            instr = instrumentList{iInst};
            idx   = indexArrays.(instr)(f_idx);

            thisFileName = allInstrumentData.(instr).Files{idx};
            thisFullPath = fullfile(rootBaseDir, upper(instr), 'Data', thisFileName);
            varName      = allInstrumentData.(instr).VarName;

            % Output folder for single-wavelet results for that instrument
            singleOutDir = fullfile(rootBaseDir, upper(instr), ...
                'Wavelet_Results', sprintf('Frame_%s', frameDateStr));
            if ~exist(singleOutDir, 'dir')
                mkdir(singleOutDir);
            end

            % Check if the single-wavelet netCDF already exists
            singleNcFile = fullfile(singleOutDir, sprintf('FrameWavelet_%s.nc', frameDateStr));
            if isfile(singleNcFile)
                % If it already exists, we could load it if we need the spec_full, 
                % but presumably we do NOT want to store huge spec_full in netCDF. 
                % For cross-wavelet, we need the actual wavelet coefficients. 
                % => We might either re-compute or store them in memory anyway.
                fprintf('Single-wavelet NetCDF for %s at %s already exists. Will re-calc wavelet in memory.\n', ...
                    instr, frameDateStr);
            end
    
            % Read data & preprocess
            data = double(ncread(thisFullPath, varName));
            data_pre = preprocessFrame(data, instr, getPreprocessingMethod(instr));
            data_filt = data_pre; %Kept in memory for the annotation function

            % Resizing
            if shrinkfactor ~= 1
                data_pre = imresize(data_pre, invshrinkfactor);
            end

            % (Optional) Radial window
            if doRadialWindow
                data_pre = applyRadialWindow(data_pre, radius_factor, decay_rate);
            end

            [rowsF, colsF] = size(data_pre);

            % Build squares if not already done
            % We can do it just for the first instrument, but let's do it for the largest dimension, etc.
            if (iInst == 1)
                x_buffer_range = (window_buffer+1) : (colsF - window_buffer);
                y_buffer_range = (window_buffer+1) : (rowsF - window_buffer);
                adjusted_frame_width  = length(x_buffer_range);
                adjusted_frame_height = length(y_buffer_range);

                pxPerSquare = round(square_size_deg / effective_deg_per_px);
                num_squares_x = ceil(adjusted_frame_width  / pxPerSquare);
                num_squares_y = ceil(adjusted_frame_height / pxPerSquare);

                squares = [];
                idxS = 1;
                for iy = 1:num_squares_y
                    for ix = 1:num_squares_x
                        x_start = floor((ix - 1) * adjusted_frame_width  / num_squares_x) + 1;
                        y_start = floor((iy - 1) * adjusted_frame_height / num_squares_y) + 1;
                        x_end   = floor(ix * adjusted_frame_width  / num_squares_x);
                        y_end   = floor(iy * adjusted_frame_height / num_squares_y);

                        if x_end > x_start && y_end > y_start
                            squares(idxS).x_range = x_buffer_range(x_start:x_end);
                            squares(idxS).y_range = y_buffer_range(y_start:y_end);
                            squares(idxS).index   = idxS;
                            idxS = idxS + 1;
                        end
                    end
                end
                squaresInfo = squares;
            else
                % Optionally, check if sizes match across instruments
                % If IR and VIS differ in dimension, you'd have to handle it here
            end

            % Compute 2D wavelet
            waveStruct = cwtft2(data_pre, 'wavelet','cauchy','scales',Scales,'angles',Angles);
            spec_full  = squeeze(waveStruct.cfs);  % [Ny, Nx, NSCALES, NANGLES]

            % Wavelet normalization
            for iS = 1:NSCALES
                spec_full(:,:,iS,:) = spec_full(:,:,iS,:) * (2 / Scales(iS));
            end

            % Average Wavelet over Scales and Angles
            % Initialize a new spec_full_avg with dimensions based on the number of squares
            nScales = size(spec_full, 3);  % Number of scales
            nAngles = size(spec_full, 4); % Number of angles
            
            spec_full_avg = zeros(num_squares_y,num_squares_x, nScales, nAngles); % (NbsquareCol*NbsquareLine x Scales x Angles)
            
            for iy = 1:num_squares_y
                for ix = 1:num_squares_x
                    % Determine the ranges for the current square
                    x_start = floor((ix - 1) * adjusted_frame_width / num_squares_x) + 1;
                    x_end = floor(ix * adjusted_frame_width / num_squares_x);
                    y_start = floor((iy - 1) * adjusted_frame_height / num_squares_y) + 1;
                    y_end = floor(iy * adjusted_frame_height / num_squares_y);
            
                    % Extract the region of spec_full corresponding to this square
                    square_region = spec_full(y_buffer_range(y_start:y_end), ...
                                               x_buffer_range(x_start:x_end), :, :);
                    
                    % Compute the average over this square for each scale and angle
                    spec_full_avg(iy, ix, :, :) = squeeze(mean(mean(square_region, 1, 'omitnan'), 2, 'omitnan'));
                end
            end

            % Save to NetCDF if not already present
            if ~isfile(singleNcFile)
                % Create
                nccreate(singleNcFile, 'spec_full_avg', ...
                    'Dimensions', {'squares_y', num_squares_y, ...
                                   'squares_x', num_squares_x, ...
                                   'scales', NSCALES, ...
                                   'angles', NANGLES}, ...
                    'Datatype', 'double');
                nccreate(singleNcFile, 'scales', 'Dimensions', {'scales', NSCALES}, 'Datatype', 'double');
                nccreate(singleNcFile, 'angles', 'Dimensions', {'angles', NANGLES}, 'Datatype', 'double');
                
                % Write
                ncwrite(singleNcFile, 'spec_full_avg', spec_full_avg);
                ncwrite(singleNcFile, 'scales', Scales);
                ncwrite(singleNcFile, 'angles', Angles);

                fprintf('Saved single-wavelet NetCDF for %s at %s -> %s\n',...
                    instr, frameDateStr, singleNcFile);

                produceAnnotatedImages(instr,spec_full, data_filt, squares, Scales, Angles, singleOutDir, ...
                    frameDateStr, pixel_size_km, 1, 1);
            else
                fprintf('NetCDF file for %s at %s already exists, skipping creation.\n', instr, frameDateStr);
            end

            % Save results in memory to do cross-wavelet later
            waveletSpecs.(instr) = spec_full;

        end % loop over each instrument

        % Step 2: Compute cross-wavelets pairwise among all instruments
        %------------------------------------------------------------------
        if numel(instrumentList) < 2
            continue;  % no cross-wavelet if only 1 instrument
        end

        % Build all pairs in alphabetical (canonical) order
        allPairs = nchoosek(sort(instrumentList), 2); 
        % allPairs is a cell array Nx2 if instrumentList is a cell array of strings

        for p = 1:size(allPairs,1)
            instA = allPairs{p,1};
            instB = allPairs{p,2};

            % Sort again to be absolutely sure we have canonical naming
            sortedPair = sort({instA, instB});
            pairName   = upper(strjoin(sortedPair, '_'));  % e.g. 'IR_VIS'

            % Create an output directory for cross-wavelet
            crossOutDir = fullfile(rootBaseDir, pairName, ...
                'CrossWavelet_Results', sprintf('Frame_%s', frameDateStr));
            if ~exist(crossOutDir, 'dir')
                mkdir(crossOutDir);
            end

            % Build canonical NetCDF name
            crossNcName = sprintf('CrossWavelet_%s.nc', frameDateStr);
            crossNcFile = fullfile(crossOutDir, crossNcName);

            if isfile(crossNcFile)
                fprintf('Cross-wavelet netCDF already exists for pair %s at %s -> skip.\n',...
                    pairName, frameDateStr);
                continue;
            end

            % Retrieve wavelet specs from memory
            specA = waveletSpecs.(instA);  % [Ny, Nx, nScales, nAngles]
            specB = waveletSpecs.(instB);

            % Compute cross
            crossSpec = specA .* conj(specB);
            % You might want to compute coherence = abs(crossSpec), phase = angle(crossSpec), etc.

            % Average over squares
            [~, ~, nScales, nAngles] = size(specA); 
            crossSpec_avg = zeros(num_squares_y, num_squares_x, nScales, nAngles);

            for iy = 1:num_squares_y
                for ix = 1:num_squares_x
                    % Determine the ranges for the current square
                    x_start = floor((ix - 1) * adjusted_frame_width / num_squares_x) + 1;
                    x_end = floor(ix * adjusted_frame_width / num_squares_x);
                    y_start = floor((iy - 1) * adjusted_frame_height / num_squares_y) + 1;
                    y_end = floor(iy * adjusted_frame_height / num_squares_y);
            
                    % Extract the region of spec_full corresponding to this square
                    square_region = spec_full(y_buffer_range(y_start:y_end), ...
                                               x_buffer_range(x_start:x_end), :, :);
                    
                    % Compute the average over this square for each scale and angle
                    crossSpec_avg(iy, ix, :, :) = squeeze(mean(mean(square_region, 1, 'omitnan'), 2, 'omitnan'));
                end
            end

            % Save cross-wavelet to NetCDF
            nccreate(crossNcFile, 'crossSpec_avg', ...
                'Dimensions', {'squares_y', num_squares_y, ...
                               'squares_x', num_squares_x, ...
                               'scales', nScales, ...
                               'angles', nAngles}, ...
                'Datatype', 'double');
            nccreate(crossNcFile, 'scales', 'Dimensions', {'scales', nScales}, 'Datatype', 'double');
            nccreate(crossNcFile, 'angles', 'Dimensions', {'angles', nAngles}, 'Datatype', 'double');

            ncwrite(crossNcFile, 'crossSpec_avg', crossSpec_avg);
            ncwrite(crossNcFile, 'scales', Scales);
            ncwrite(crossNcFile, 'angles', Angles);

            fprintf('Saved cross-wavelet for pair %s at %s -> %s\n', ...
                pairName, frameDateStr, crossNcFile);

            % Optionally produce images for cross-spectra (phase, coherence, etc.)
            % e.g. produceAnnotatedImages(...), specialized for crossSpec
        end
    end

    fprintf('\nAll done. Single & cross-wavelet computations complete.\n');
end
%% Helper Functions
%==============================================================================
function [fileNames, fileTimestamps, variableName] = getDateRangeFiles(dataDir, startDate, endDate)
    % Suppose que les fichiers sont du type:
    % IR_YYYY_MM_DD_HH_MM.nc

    ncFiles = dir(fullfile(dataDir, '*.nc'));
    if isempty(ncFiles)
        fileNames = {};
        fileTimestamps = [];
        variableName = '';
        warning('No .nc files in %s', dataDir);
        return;
    end

    % Pour la variable, on check dans le 1er fichier
    firstFile = fullfile(ncFiles(1).folder, ncFiles(1).name);
    info = ncinfo(firstFile);
    varList = {info.Variables.Name};
    if any(strcmp(varList, 'CMI'))
        variableName = 'CMI';
    elseif any(strcmp(varList, 'Rad'))
        variableName = 'Rad';
    else
        variableName = varList{1};
        warning('No standard var found; using %s', variableName);
    end

    fileNames = {};
    fileTimestamps = datetime([], 'ConvertFrom', 'datenum'); % Ensure it is initialized as an empty datetime array

    for i = 1:numel(ncFiles)
        fn = ncFiles(i).name;
        % fn = e.g. IR_2023_10_12_01_15.nc
        parts = split(fn, '_');
        if numel(parts) < 6
            continue;
        end
        % parts{1} = 'IR'
        % parts{2} = '2023'
        % parts{3} = '10'
        % parts{4} = '12'
        % parts{5} = '01'
        % parts{6} = '15.nc'

        yyyy = str2double(parts{2});
        mm   = str2double(parts{3});
        dd   = str2double(parts{4});
        HH   = str2double(parts{5});
        minPart = erase(parts{6}, '.nc'); 
        MN   = str2double(minPart);

        try
            thisDT = datetime(yyyy, mm, dd, HH, MN, 0);
        catch
            continue;
        end

        if thisDT >= startDate && thisDT <= endDate
            fileNames{end+1} = fn;
            fileTimestamps(end+1) = thisDT;
        end
    end

    [fileTimestamps, idxSort] = sort(fileTimestamps);
    fileNames = fileNames(idxSort);
end

function data_preprocessed = preprocessFrame(data, dataType, methodName)
    % dataType: e.g., 'IR' or 'VIS'
    % methodName: e.g., 'highpass_50_sqrt', 'raw_normalized', etc.

    switch upper(dataType)
        case 'IR'
            %-----------------------------------------------------
            % (1) IR-specific threshold or masking
            %-----------------------------------------------------
            threshold = prctile(data(:), 5);
            data(data < threshold) = NaN; 
            
            % Remember which pixels were set to NaN
            nan_mask = isnan(data);
            
            % Fill those NaNs with some fallback (mean of valid data)
            data(nan_mask) = mean(data,'all','omitnan');

            %-----------------------------------------------------
            % (2) Next, call your method-based process
            %-----------------------------------------------------
            data_preprocessed = processDataMethod(data, methodName);

            %-----------------------------------------------------
            % (3) If you want to “clamp back” to some percentile
            %     for those NaN locations, do that now:
            %-----------------------------------------------------
            % For instance, if you previously did:
            %   data_pro(nan_mask') = prctile(data_pro(:), 50);
            % replicate it here:
            fill_value = prctile(data_preprocessed(:), 50);
            data_preprocessed(nan_mask') = fill_value;

        case 'VIS'
            %-----------------------------------------------------
            % (1) VIS-specific threshold or dynamic range
            %-----------------------------------------------------
            lowerBound = prctile(data(:), 10);
            upperBound = prctile(data(:), 99);
            data(data < lowerBound) = lowerBound;
            data(data > upperBound) = upperBound;

            % (no special NaN masking for VIS?)

            %-----------------------------------------------------
            % (2) Next, call your method-based process
            %-----------------------------------------------------
            data_preprocessed = processDataMethod(data, methodName);

        otherwise
            %-----------------------------------------------------
            % If unrecognized, you can either throw an error or do
            % a fallback approach
            %-----------------------------------------------------
            warning('Unrecognized dataType: %s. Using fallback.', dataType);
            data_preprocessed = processDataMethod(data, methodName);
    end
end

function img_processed = processDataMethod(data, methodName)
    switch methodName
        case 'none'
            img_processed=data';
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

function img_out = applyHighPass(data, filterWidth, doSqrtEnhance)
    lowPass = imgaussfilt(data, filterWidth);
    highPass = data - lowPass;
    if doSqrtEnhance
        highPass = sqrt(abs(highPass)) .* sign(highPass);
    end
    clipMin = -3;
    clipMax = 3;
    highPass(highPass < clipMin) = clipMin;
    highPass(highPass > clipMax) = clipMax;
    img_out = (highPass - clipMin) / (clipMax - clipMin);
    img_out = 1 - img_out;
end

function data_win = applyRadialWindow(data_in, radius_factor, decay_rate)
    [rows, cols] = size(data_in);
    cx = cols/2; 
    cy = rows/2;
    [X, Y] = meshgrid(1:cols, 1:rows);
    R = sqrt((X - cx).^2 + (Y - cy).^2);

    maxR  = radius_factor * min(cx, cy);
    window = 1 ./ (1 + exp(decay_rate * (R - maxR)));

    data_win = data_in .* window;
end

function produceAnnotatedImages(dataType, spec_full, data_background, squares, ...
                                Scales, Angles, outDir, frameDateStr, ...
                                pixel_size_km, clevfactor, saverose)
% PRODUCEANNOTATEDIMAGES Main processing function for wavelet analysis visualization.
%
% Inputs:
%   dataType           - String type: 'IR' or 'VIS'
%   spec_full          - 4D wavelet coefficients [Ny_sh, Nx_sh, nScales, nAngles]
%   data_background    - Original 2D background data
%   squares            - Structure array of ROI information
%   Scales             - Vector of wavelet scales
%   Angles             - Vector of wavelet angles
%   outDir             - Output directory path
%   frameDateStr       - Frame timestamp string
%   pixel_size_km      - Pixel size in kilometers
%   clevfactor         - Contour level adjustment factor
%   saverose           - Boolean flag to save (true) or skip saving (false) the waverose image
%
% -------------------------------------------------------------------------
% 1) Dimensions & Initial Processing
% -------------------------------------------------------------------------
[Ny_sh, Nx_sh, nScales, nAngles] = size(spec_full);
[Ny_orig, Nx_orig] = size(data_background);

% Compute wavelet power: power = |coefficients|^2
power = abs(spec_full).^2;

% Compute average power over spatial dimensions (innerpower is scale-angle power)
innerpower = squeeze(mean(mean(power, 1, 'omitnan'), 2, 'omitnan'));

% -------------------------------------------------------------------------
% 2) Wave-Rose Detection and Visualization
% -------------------------------------------------------------------------
% 2.1) Interpolate power to finer grid for smoother rose
nAngles_fine = 4 * nAngles;
nScales_fine = 4 * nScales;
Angles_fine = linspace(min(Angles), max(Angles), nAngles_fine);
Scales_fine = linspace(min(Scales), max(Scales), nScales_fine);

% Create mesh for original grid
[Theta_orig, R_orig] = meshgrid(Angles, Scales);

% Create mesh for finer grid
[Theta_fine, R_fine] = meshgrid(Angles_fine, Scales_fine);

% Interpolate using a gridded interpolant
F = griddedInterpolant(Theta_orig', R_orig', innerpower', 'spline');
innerpower_fine = F(Theta_fine', R_fine')';

% Convert polar to Cartesian for positive hemisphere (ax1)
[X_pos_fine, Y_pos_fine] = pol2cart(Theta_fine, R_fine);

% Convert polar to Cartesian for negative hemisphere (ax2) by adding pi
[X_neg_fine, Y_neg_fine] = pol2cart(Theta_fine + pi, R_fine);

% Create figure and two axes (for dual-hemisphere plotting)
figRose = figure('visible','off');
ax1 = axes('Position',[0.1 0.1 0.75 0.75]); 
hold(ax1, 'on');
pcolor(ax1, X_pos_fine, Y_pos_fine, innerpower_fine);
shading(ax1, 'interp');
colormap(ax1, 'parula');
axis(ax1, 'equal', 'tight', 'off');

ax2 = axes('Position', ax1.Position, 'Color','none', 'HitTest','off');
hold(ax2, 'on');
pcolor(ax2, X_neg_fine, Y_neg_fine, innerpower_fine);
shading(ax2, 'interp');
axis(ax2, 'equal', 'tight', 'off');

% Ensure ax1 stays on top for interactions, link properties
uistack(ax1, 'top');
linkprop([ax1 ax2], {'XLim','YLim','Position','CameraPosition','CameraUpVector'});

% 2.2) Threshold-based peak detection on original (non-interpolated) data
threshold_orig = 0.7 * max(innerpower(:));
bwMask_orig = innerpower >= threshold_orig;
CC = bwconncomp(bwMask_orig, 4);
numPeaks = CC.NumObjects;

% Contour on interpolated data for visualization
threshold_viz = 0.7 * max(innerpower_fine(:));
contour(ax1, X_pos_fine, Y_pos_fine, innerpower_fine, [threshold_viz threshold_viz], 'r-', 'LineWidth',2);
contour(ax2, X_neg_fine, Y_neg_fine, innerpower_fine, [threshold_viz threshold_viz], 'r-', 'LineWidth',2);

%% 2.2) Label peaks in the rose
% We compute each region’s mean scale and angle, then place a text label.
for pk = 1:numPeaks
    % Convert linear indices to subscripts (scaleIndex, angleIndex)
    [scaleIndices, angleIndices] = ind2sub(size(bwMask_orig), CC.PixelIdxList{pk});
    
    % Compute mean scale and mean angle for that region
    meanScale = mean(Scales(scaleIndices), 'omitnan');
    meanAngle = mean(Angles(angleIndices), 'omitnan');
    
    % Depending on the meanAngle, place the label in ax1 or ax2
    [x_peak, y_peak] = pol2cart(meanAngle, meanScale);
    text(ax1, x_peak, y_peak, sprintf('%d', pk), ...
        'Color','k', 'FontWeight','bold', ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
end

% -------------------------------------------------------------------------
% 3) Annotations and Labels
% -------------------------------------------------------------------------
% 3.1) Draw radial grid circles for scales
if exist('pixel_size_km', 'var')
    wavelengths_km = Scales * pi/sqrt(2) * pixel_size_km;
else
    wavelengths_km = Scales;
end

for i = 1:length(Scales)
    theta_ring = linspace(0, 2*pi, 100);
    [x_ring, y_ring] = pol2cart(theta_ring, Scales(i));
    plot(ax1, x_ring, y_ring, 'k--', 'LineWidth', 0.5);
    plot(ax2, x_ring, y_ring, 'k--', 'LineWidth', 0.5);
    text(ax1, Scales(i)*1.05, 0, sprintf('%.1f km', wavelengths_km(i)),...
         'HorizontalAlignment', 'left', 'FontSize', 8);
end

% 3.2) Draw angular lines and labels
angle_ticks = linspace(0, 2*pi, 13);
angle_labels = {'0','\pi/6','\pi/3','\pi/2','2\pi/3','5\pi/6','\pi',...
                '7\pi/6','4\pi/3','3\pi/2','5\pi/3','11\pi/6','2\pi'};
max_r = max(Scales)*1.1;

for i = 1:length(angle_ticks)
    [x_label, y_label] = pol2cart(angle_ticks(i), max_r);
    line(ax1, [0 x_label], [0 y_label], 'Color',[0.5 0.5 0.5], 'LineStyle','--');
    line(ax2, [0 x_label], [0 y_label], 'Color',[0.5 0.5 0.5], 'LineStyle','--');
    
    if angle_ticks(i) <= pi
        text(ax1, x_label*1.05, y_label*1.05, angle_labels{i},...
             'HorizontalAlignment','center', 'FontSize',8);
    else
        text(ax2, x_label*1.05, y_label*1.05, angle_labels{i},...
             'HorizontalAlignment','center', 'FontSize',8);
    end
end

% 3.3) Colorbar
c = colorbar(ax1, 'Location', 'eastoutside');
c.Label.String = 'Wavelet Power';
c.Label.FontWeight = 'bold';
ax1_pos = ax1.Position;
ax1_pos(3) = ax1_pos(3) * 0.85;  % Make room for colorbar
ax1.Position = ax1_pos;
ax2.Position = ax1_pos;

% 3.4) Add final title
title(ax1, sprintf('Polar Wave-Rose: %d Significant Regions', numPeaks),...
      'FontSize',12, 'FontWeight','bold');

% 3.5) Conditional saving of the wave-rose
if saverose
    roseName = fullfile(outDir, sprintf('WaveRose_%s.png', frameDateStr));
    exportgraphics(figRose, roseName, 'Resolution', 300);
end
close(figRose);

% -------------------------------------------------------------------------
% 4) Region Processing and Image Output
% -------------------------------------------------------------------------
% Compute scale factor between original and wavelet grid
scaleFactorX = Nx_orig / Nx_sh;
scaleFactorY = Ny_orig / Ny_sh;

% Prepare a cell array storing scale/angle info for each region
peakRegions = cell(numPeaks,1);
for pk = 1:numPeaks
    [scaleIndices, angleIndices] = ind2sub(size(bwMask_orig), CC.PixelIdxList{pk});
    
    % Convert scale index to kilometers (if pixel_size_km is provided)
    scales_km = Scales(scaleIndices) * pi/sqrt(2) * pixel_size_km;
    angles_deg = rad2deg(Angles(angleIndices));
    
    % Format as strings
    scale_str = join(split(num2str(scales_km,'%.1f ')), '/');
    angle_str = join(split(num2str(angles_deg,'%.0f ')), '/');
    
    peakRegions{pk} = struct(...
        'ScaleIndices', scaleIndices, ...
        'AngleIndices', angleIndices, ...
        'ScaleStr', scale_str{1}, ...
        'AngleStr', angle_str{1});
end

% Generate output images for each peak
for pk = 1:numPeaks
    % Sum real parts of coefficients in the current region, and also accumulate power
    waveSum = zeros(Ny_sh, Nx_sh);
    wavePower = zeros(Ny_sh, Nx_sh);
    currentRegion = peakRegions{pk};
    
    for jj = 1:numel(currentRegion.ScaleIndices)
        s_idx = currentRegion.ScaleIndices(jj);
        a_idx = currentRegion.AngleIndices(jj);
        coeff = spec_full(:,:,s_idx,a_idx);
        waveSum = waveSum + real(coeff);
        wavePower = wavePower + abs(coeff).^2;
    end
    
    % Upsample to original image size
    waveSum_up = imresize(waveSum, [Ny_orig, Nx_orig]);
    wavePower_up = imresize(wavePower, [Ny_orig, Nx_orig]);
    
    % Create a figure for the overlay
    fig = figure('visible','off');
    switch upper(dataType)
        case 'IR'
            imagesc(data_background, [0 1])
        case 'VIS'
            image(data_background);
        otherwise
            error('Unknown dataType.');
    end
    colormap(gray);
    axis image off;
    hold on;
    
    % Define contour levels
    num_levels = 5; % example number of levels
    std_real = std(waveSum_up(:));
    pos_levels = linspace(std_real, max(waveSum_up(:)), num_levels) / clevfactor;
    neg_levels = linspace(-std_real, min(waveSum_up(:)), num_levels) / clevfactor;
    
    std_power = std(wavePower_up(:));
    power_levels = linspace(std_power, max(wavePower_up(:)), num_levels) / clevfactor^2;
    
    % Draw contours for waveSum (positive in red, negative in blue) and wavePower (white)
    contour(waveSum_up, pos_levels, 'LineColor','red',  'LineWidth',0.5);
    contour(waveSum_up, neg_levels, 'LineColor','blue', 'LineWidth',0.5);
    contour(wavePower_up, power_levels, 'LineColor',[0.99 0.99 0.99], 'LineWidth',0.5);
    
    % Draw the ROI rectangles
    for sq = 1:numel(squares)
        xPos_orig = squares(sq).x_range(1) * scaleFactorX;
        yPos_orig = squares(sq).y_range(1) * scaleFactorY;
        w_orig = length(squares(sq).x_range) * scaleFactorX;
        h_orig = length(squares(sq).y_range) * scaleFactorY;
        rectangle('Position',[xPos_orig, yPos_orig, w_orig, h_orig],...
                  'EdgeColor','k','LineWidth',1);
    end
    
    % Multi-line title with scales and angles
    titleText = {sprintf('Instrument X - %s', frameDateStr), ...
                 sprintf('Peak %d/%d - Scales: %s km', pk, numPeaks, peakRegions{pk}.ScaleStr), ...
                 sprintf('Angles: %s°', peakRegions{pk}.AngleStr)};
    title(titleText, 'Color', 'k', 'FontWeight','bold', 'FontSize', 10, 'Interpreter', 'none');
    
    % Save and close
    outName = fullfile(outDir, sprintf('Frame_%s_Region%02d.png', frameDateStr, pk));
    saveas(fig, outName);
    close(fig);
end

end

function renameAndOrganizeFiles(dataType, startDate, endDate, sourceRootDir, outRootDir)
% RENAMEANDORGANIZEFILES
% 
% Scans a source directory, searches for .nc files,
% extracts their timestamps, checks whether they lie in the [startDate, endDate] range,
% rounds that timestamp to the nearest quarter-hour, then renames each file to
% "INSTRUMENT_YYYY_MM_DD_HH_MM.nc" (using the rounded time).
% Finally, it copies these renamed files into:
%   outRootDir\<dataType>\Data
%
% INPUTS:
%   dataType      = 'IR' or 'VIS'
%   startDate     = datetime (begin)
%   endDate       = datetime (end)
%   sourceRootDir = e.g. 'C:\Users\admin\Box\GOES2go_satellite_downloads'
%                   (we expect subfolders like '4km_SEPAC_IR', etc.)
%   outRootDir    = e.g. 'C:\Users\admin\Box\GWaves_2023_10_11-14_SEPAC'
%
% This function copies files instead of moving them, so the originals remain in place.

    % 1) Determine the specific source folder based on the instrument
    switch upper(dataType)
        case 'IR'
            originalDir = fullfile(sourceRootDir, '4km_SEPAC_IR');
        case 'VIS'
            originalDir = fullfile(sourceRootDir, '4km_SEPAC_VIS');
        otherwise
            error('Unknown dataType: %s', dataType);
    end

    if ~exist(originalDir, 'dir')
        error('Source directory does not exist: %s', originalDir);
    end

    % 2) Create the destination folder under outRootDir\<dataType>\Data
    destDir = fullfile(outRootDir, upper(dataType), 'Data');
    if ~exist(destDir, 'dir')
        mkdir(destDir);
    end

    % 3) List all .nc files in the source directory
    ncFiles = dir(fullfile(originalDir, '*.nc*'));
    if isempty(ncFiles)
        warning('No .nc files found in %s', originalDir);
        return;
    end

    % 4) Loop through all found files
    for i = 1:numel(ncFiles)
        oldName = ncFiles(i).name;  % e.g. OR_ABI-L2-CMIPC-M6C13_G16_s2023245...
        oldPath = fullfile(originalDir, oldName);

        % Extract the timestamp between "_s" and "_e"
        tsStr = extractBetween(oldName, '_s', '_e');
        if isempty(tsStr)
            continue;
        end
        tsStr = tsStr{1};

        % Convert the extracted string to a datetime
        try
            fileTS = datetime(tsStr, 'InputFormat', 'uuuuDDDHHmmssSSS');
        catch
            continue; % unrecognized format
        end

        % Check if the timestamp lies within [startDate, endDate]
        if fileTS < startDate || fileTS > endDate
            continue;
        end

        % 5) Round timestamp to the nearest quarter hour
        fileTS_rounded = roundToQuarterHour(fileTS);

        % 6) Construct the new name: INSTRUMENT_YYYY_MM_DD_HH_MM.nc
        YYYY = year(fileTS_rounded);
        MM   = month(fileTS_rounded);
        DD   = day(fileTS_rounded);
        HH   = hour(fileTS_rounded);
        MN   = minute(fileTS_rounded);

        newName = sprintf('%s_%04d_%02d_%02d_%02d_%02d.nc', ...
                          upper(dataType), YYYY, MM, DD, HH, MN);

        newPath = fullfile(destDir, newName);

        % 7) Copy the file under the new name
        if ~isfile(newPath)
            copyfile(oldPath, newPath);
            fprintf('Copied: %s -> %s\n', oldName, newName);
        else
            % If the file already exists, you may overwrite or skip
            fprintf('[!] File %s already exists. Skipping.\n', newName);
        end
    end
end

%% Subfunction: rounding to the nearest quarter hour
function roundedDT = roundToQuarterHour(originalDT)
    % roundToQuarterHour
    % Rounds a datetime to the nearest 15-minute interval.
    % Example: 12:07 -> 12:00, 12:08 -> 12:15

    minutesFromQuarter = mod(minute(originalDT), 15);
    if minutesFromQuarter < 7.5
        % closer to the previous quarter hour
        roundedDT = dateshift(originalDT, 'start', 'minute') - minutes(minutesFromQuarter);
    else
        % closer to the next quarter hour
        roundedDT = dateshift(originalDT, 'start', 'minute') + minutes(15 - minutesFromQuarter);
    end
end

function methodName = getPreprocessingMethod(dataType)
    switch upper(dataType)
        case 'IR'
            methodName = 'highpass_50_sqrt';
        case 'VIS'
            methodName = 'none';
        otherwise
            error('Unknown dataType.');
    end
end


