%%

% Exemple:
%startDate = datetime(2023,09,28,1,0,0);
%endDate   = datetime(2023,09,29,1,0,0);
startDate = datetime(2023, 10, 11, 1, 0, 0);
endDate   = datetime(2023, 10, 14, 9, 0, 0); 


sourceRoot   = 'C:\Users\admin\Box\GOES2go_satellite_downloads';
outputRoot   = 'C:\Users\admin\Box\GWaves_2023_10_11-14_SEPAC';

% For IR
renameAndOrganizeFiles('IR', startDate, endDate, sourceRoot, outputRoot);

% For VIS
renameAndOrganizeFiles('VIS', startDate, endDate, sourceRoot, outputRoot);

SingleFrameWavelet_Explicit_

function SingleFrameWavelet_Explicit_
%==========================================================================
% SINGLEFRAMEWAVELET_CODE
%
% This function processes GOES satellite data (IR or VIS) in a single-frame
% context. It loads .nc files within a specified date/time range, performs
% preprocessing, wavelet transforms, region-based averaging, and generates
% annotated images/plots of gravity-wave-like signals.
%
% All parameters and thresholds are defined in the "Variables Setup" section.
%==========================================================================
%% 1) VARIABLES SETUP
%--------------------------------------------------------------------------

% ----------- DATE/TIME & DATA-TYPE SETTINGS ------------------------------
startDate              = datetime(2023, 10, 11, 1, 0, 0);   % [Used in getDateRangeFiles] Start of date range
endDate                = datetime(2023, 10, 14, 9, 0, 0);   % [Used in getDateRangeFiles] End of date range
dataType               = 'IR';                             % [Used throughout code] 'IR' or 'VIS' determines paths & thresholds


% ----------- FOLDER/PATH SETTINGS ----------------------------------------
rootSepacDir           = 'C:\Users\admin\Box\GWaves_2023_10_11-14_SEPAC';   % [Used to build waveletResultsDir, etc.]
boxRootDir_IR          = fullfile(rootSepacDir, 'IR', 'Data');             % [Used if dataType='IR'] Path to IR .nc files
boxRootDir_VIS         = fullfile(rootSepacDir, 'VIS', 'Data');            % [Used if dataType='VIS'] Path to VIS .nc files
waveletResultsDir      = '';                                               % [Will be set below based on dataType]

% ----------- SPATIAL SCALING & RESIZING ----------------------------------
degrees_per_pixel      = 0.04;   % [Used in building scale. 0.04° ~ typical GOES 4km scale]
km_per_degree          = 111.32; % [Used in building scale. 1 degree ~ 111.32 km at Earth surface approx.]
shrinkfactor           = 2;      % [Used in imresize step. 2 => final image has half the original resolution]
invshrinkfactor        = 1/shrinkfactor;  % [Used in imresize for convenience]
original_px_km         = degrees_per_pixel * km_per_degree;  % [Intermediate: original ~4.4528 km/px if 0.04° used]
pixel_size_km          = original_px_km * shrinkfactor;       % [Resulting px size after shrink]

% ----------- WAVELET PARAMETERS ------------------------------------------
Angles                 = pi/9 : pi/9 : (8*pi/9); % [Used in cwtft2. This sets wavelet angle sampling ~20° increments]
Scales                 = [2, 4, 8, 16, 32, 64, 128]; % [Used in cwtft2. This sets wavelet scale in "pixel" units]
NANGLES                = numel(Angles);   % [For reference only]
NSCALES                = numel(Scales);   % [For reference only]

% ----------- WINDOWING (OPTIONAL) ----------------------------------------
doWindow               = true;    % [Used in main loop to apply a window function if true]
windowType             = 'rectangular';  % ['rectangular' or 'radial', used in applyRectangularWindow/applyRadialWindow]
radius_factor          = 0.8;     % [Used inside window function to define where window decays]
decay_rate             = 10;      % [Used inside window function controlling how steep the window edge is]

% ----------- SQUARE-PARTITIONING PARAMETERS ------------------------------
window_buffer          = 10;   % [Used to avoid edges. We keep an inner region by skipping 'window_buffer' px at each edge]
square_size_deg        = 5;    % [Used to define about 5° squares. This is converted to px for partitioning]

% ----------- IR PREPROCESSING THRESHOLDS ---------------------------------
IR_threshold           = 280;  % [Used in preprocessFrame for IR. Values < 275 => set to NaN -> then filled]
IR_fillPercentile      = 50;   % [Used after IR_threshold masking. We fill masked areas with this percentile of final data]

% ----------- VIS PREPROCESSING THRESHOLDS --------------------------------
VIS_lowerPercentile    = 10;   % [Used in preprocessFrame for VIS. Lower bound clip at prctile(..., 10)]
VIS_upperPercentile    = 99;   % [Used in preprocessFrame for VIS. Upper bound clip at prctile(..., 99)]
VIS_fillPercentile     = 50;    % Fill VIS NaN pixels with this percentile

% ----------- HIGH-PASS FILTER SETTINGS -----------------------------------
clipMinHP              = -3;    % [Used in applyHighPass to clamp negative extremes of highpass signal]
clipMaxHP              =  3;    % [Used in applyHighPass to clamp positive extremes of highpass signal]
lowPassFilterWidth_20  = 20;    % [Used in applyHighPass('highpass_20'...)]
lowPassFilterWidth_50  = 50;    % [Used in applyHighPass('highpass_50'...)]
lowPassFilterWidth_100 = 100;   % [Used in applyHighPass('highpass_100'...)]

% ----------- ALTERNATE PREPROCESSING MODES -------------------------------
IR_methodName          = 'highpass_50_sqrt';  % [Default method if dataType='IR']
VIS_methodName         = 'none';              % [Default method if dataType='VIS']

% ----------- WAVE-ROSE & PEAK DETECTION ----------------------------------
nAngles_fineFactor     = 4;            % [Used in produceAnnotatedImages to refine angular resolution for rose plot]
nScales_fineFactor     = 4;            % [Used in produceAnnotatedImages to refine radial resolution for rose plot]
peakDetectionFactor    = 0.1;          % [Used for threshold_orig = mean(...) + 0.1*std(...)]
contourArray           = [95 97 99];   %[Used as either percentiles or absolute values for the contouring]
ArrayMode              = 'percentile'; % determines if the array represent 'percentile' or 'absolute' values

% ----------- IMAGE ANNOTATIONS & OUTPUT ----------------------------------
saverose               = 1;  % [Used in produceAnnotatedImages: if 1 => saves the wave-rose figure to disk]

%==========================================================================

%% 2) BASIC SETUP (SELECT DIRECTORIES BASED ON DATATYPE)
switch upper(dataType)
    case 'IR'
        boxRootDir = boxRootDir_IR;
        waveletResultsDir = fullfile(rootSepacDir, 'IR', 'Wavelet_Results');
        methodName = IR_methodName; 
    case 'VIS'
        boxRootDir = boxRootDir_VIS;
        waveletResultsDir = fullfile(rootSepacDir, 'VIS', 'Wavelet_Results');
        methodName = VIS_methodName;
    otherwise
        error('Unknown dataType.');
end

% Check existence / create wavelet results dir
if ~exist(boxRootDir, 'dir')
    error('Box root directory does not exist: %s', boxRootDir);
end
if ~exist(waveletResultsDir, 'dir')
    mkdir(waveletResultsDir);
end

% Retrieve all .nc in the date range
[fileNames, fileTimestamps, variableName] = getDateRangeFiles(boxRootDir, startDate, endDate);
num_frames = numel(fileNames);
if num_frames < 1
    error('No files found in the specified range.');
end
fprintf('Found %d frames from %s to %s\n', num_frames, ...
    datestr(startDate), datestr(endDate));

%% 3) LOOP OVER FRAMES (SINGLE-FRAME WAVELET + ANNOTATED IMAGES)
for f_idx = 1 : num_frames
    
    %----------------------------------------------------------------------
    % 3a) GET FRAME AND PREPROCESS
    %----------------------------------------------------------------------
    fileName = fileNames{f_idx};
    fileTime = fileTimestamps(f_idx);
    fprintf('[%d/%d] Processing file: %s\n', f_idx, num_frames, fileName);

    % Create subfolder for this frame
    frameDateStr = datestr(fileTime, 'yyyy_mm_dd_HHMMSS');
    outDir = fullfile(waveletResultsDir, sprintf('Frame_%s', frameDateStr));
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    % Read raw data
    data = double(ncread(fullfile(boxRootDir, fileName), variableName));
    fullPath = fullfile(boxRootDir, fileName);
%%
    % Preprocess (IR or VIS)
    data_pro = preprocessFrame(data, dataType, ...
        methodName,fullPath,fileTime, ...
        IR_threshold, IR_fillPercentile, ...
        VIS_lowerPercentile, VIS_upperPercentile, ...
        VIS_fillPercentile, clipMinHP, clipMaxHP, ...
        lowPassFilterWidth_20, lowPassFilterWidth_50, lowPassFilterWidth_100);

    % Keep a copy for final background overlays (no resizing)
    data_filt = data_pro;

    % Apply shrink factor if needed
    if shrinkfactor ~= 1
        data_pro = imresize(data_pro, invshrinkfactor);
    end

    % Apply optional window function
    if doWindow
        switch lower(windowType)
            case 'radial'
                data_pro = applyRadialWindow(data_pro, radius_factor, decay_rate);
            case 'rectangular'
                data_pro = applyRectangularWindow(data_pro, radius_factor, decay_rate);
            otherwise
                warning('Unknown window type: %s. No window applied.', windowType);
        end
    end
%%
    %----------------------------------------------------------------------
    % 3b) BUILD SQUARES (ROI PARTITION)
    %----------------------------------------------------------------------
    [rowsF, colsF] = size(data_pro);

    x_buffer_range = (window_buffer+1) : (colsF - window_buffer);
    y_buffer_range = (window_buffer+1) : (rowsF - window_buffer);

    adjusted_frame_width  = length(x_buffer_range);
    adjusted_frame_height = length(y_buffer_range);

    effective_deg_per_px = degrees_per_pixel; 
    square_size_px = round(square_size_deg / effective_deg_per_px);

    num_squares_x = ceil(adjusted_frame_width  / square_size_px);
    num_squares_y = ceil(adjusted_frame_height / square_size_px);

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

    %----------------------------------------------------------------------
    % 3c) WAVELET TRANSFORM + NORMALISATION
    %----------------------------------------------------------------------
    waveStruct = cwtft2(data_pro, 'wavelet','cauchy','scales',Scales,'angles',Angles);
    spec_full  = squeeze(waveStruct.cfs);  % spec_full: [Ny, Nx, nScales, nAngles]

    % Normalise wavelet coefficients by (2/Scale)
    for iS = 1:numel(Scales)
        spec_full(:,:,iS,:) = spec_full(:,:,iS,:) * (2/Scales(iS));
    end

    %----------------------------------------------------------------------
    % 3d) AVERAGE WAVELET OVER EACH ROI
    %----------------------------------------------------------------------
    [Ny_sh, Nx_sh, nScales, nAngles] = size(spec_full);
    spec_full_avg = zeros(num_squares_y,num_squares_x, nScales, nAngles);
    for iy = 1:num_squares_y
        for ix = 1:num_squares_x
            x_start = floor((ix - 1) * adjusted_frame_width / num_squares_x) + 1;
            x_end   = floor(ix * adjusted_frame_width / num_squares_x);
            y_start = floor((iy - 1) * adjusted_frame_height / num_squares_y) + 1;
            y_end   = floor(iy * adjusted_frame_height / num_squares_y);

            square_region = spec_full(y_buffer_range(y_start:y_end), ...
                                      x_buffer_range(x_start:x_end), :, :);

            spec_full_avg(iy, ix, :, :) = squeeze(mean(mean(square_region, 1, 'omitnan'), 2, 'omitnan'));
        end
    end

    %----------------------------------------------------------------------
    % 3e) SAVE RESULTS TO NETCDF (IF NOT ALREADY PRESENT)
    %----------------------------------------------------------------------
    ncFileName = fullfile(outDir, sprintf('FrameWavelet_%s.nc', frameDateStr));
    if ~isfile(ncFileName)
        nccreate(ncFileName, 'spec_full_avg', ...
            'Dimensions', {'squares_y', num_squares_y, ...
                           'squares_x', num_squares_x, ...
                           'scales', nScales, ...
                           'angles', nAngles}, ...
            'Datatype', 'double');
        nccreate(ncFileName, 'scales', 'Dimensions', {'scales', nScales}, 'Datatype', 'double');
        nccreate(ncFileName, 'angles', 'Dimensions', {'angles', nAngles}, 'Datatype', 'double');

        ncwrite(ncFileName, 'spec_full_avg', spec_full_avg);
        ncwrite(ncFileName, 'scales', Scales);
        ncwrite(ncFileName, 'angles', Angles);

        fprintf('Saved wavelet results to NetCDF: %s\n', ncFileName);
    else
        fprintf('NetCDF file already exists, skipping creation: %s\n', ncFileName);
    end

    %----------------------------------------------------------------------
    % 3f) PRODUCE ANNOTATED IMAGES (WAVE ROSE + OVERLAYS)
    %----------------------------------------------------------------------
    produceAnnotatedImages(dataType, spec_full, data_filt, squares, ...
        Scales, Angles, outDir, frameDateStr, ...
        saverose, nAngles_fineFactor, nScales_fineFactor, ...
        peakDetectionFactor,ArrayMode, contourArray)
end

fprintf('Single-frame wavelet processing completed. Results in %s\n', waveletResultsDir);
end
%==========================================================================


%% HELPER FUNCTIONS
%==========================================================================
function [fileNames, fileTimestamps, variableName] = getDateRangeFiles(dataDir, startDate, endDate)
% GETDATERANGEFILES
%  Scans dataDir for *.nc, checks timestamps in filenames of the form
%    [Type]_YYYY_MM_DD_HH_MM.nc
%  and returns arrays restricted to [startDate, endDate].

    ncFiles = dir(fullfile(dataDir, '*.nc'));
    if isempty(ncFiles)
        fileNames = {};
        fileTimestamps = [];
        variableName = '';
        warning('No .nc files in %s', dataDir);
        return;
    end

    % Identify the main data variable name (CMI, Rad, or fallback)
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
    fileTimestamps = datetime([], 'ConvertFrom', 'datenum'); 

    for i = 1:numel(ncFiles)
        fn = ncFiles(i).name;
        parts = split(fn, '_');
        if numel(parts) < 6
            continue;
        end
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
%--------------------------------------------------------------------------

function data_preprocessed = preprocessFrame(data, dataType, methodName, fullPath, ...
    fileTime, IR_threshold, IR_fillPercentile, ...
    VIS_lowerPercentile, VIS_upperPercentile, ...
    VIS_fillPercentile, clipMinHP, clipMaxHP, ...
    lpWidth20, lpWidth50, lpWidth100)
% PREPROCESSFRAME
%  Applies data-type-specific thresholds and then calls the chosen
%  method-based process. Fills or clamps as needed.

    switch upper(dataType)
        case 'IR'
            % (1) IR-specific threshold or masking
            data(data < IR_threshold) = NaN; 
            nan_mask = isnan(data);
            
            % Fill masked region with the average (or any fallback)
            data(nan_mask) = mean(data(:), 'omitnan');

            % (2) Process by methodName
            data_preprocessed = processDataMethod(data, methodName, ...
                clipMinHP, clipMaxHP, ...
                lpWidth20, lpWidth50, lpWidth100);

            % (3) Fill masked region with a chosen percentile value
            fill_value = prctile(data_preprocessed(:), IR_fillPercentile);
            data_preprocessed(nan_mask') = fill_value;

        case 'VIS'
            % (1) VIS-specific threshold or dynamic range
            lowerBound = prctile(data(:), VIS_lowerPercentile);
            upperBound = prctile(data(:), VIS_upperPercentile);
            data(data < lowerBound) = lowerBound;
            data(data > upperBound) = upperBound;

            nan_mask = isnan(data);
            
            % (2) Correct for uneven solar exposure

            lat_vec = ncread(fullPath, 'latitude');   % [1125 x 1] => vecteur lat vector
            lon_vec = ncread(fullPath, 'longitude');  % [1500 x 1] => vecteur lon vector
            data    = ncread(fullPath, 'Rad');        % [1500 x 1125] => [lon, lat]
        
            % Transpose [lat, lon]
            % data(i,j) => i=lat, j=lon
            data = data.';
            data(data < 0) = 0;
       
            [lonGrid, latGrid] = meshgrid(lon_vec, lat_vec);
            %  => size(LatGrid) = size(LonGrid) = [1125 x 1500]
        
            % We convert datetime into datenum for the function.
            dt_num = datenum(fileTime);
        
            % Solar parameters
            time_zone = 0;    
            rotation  = 0;     
            dst       = false; 

            % Getting the insolation for each pixel
            insolation_grid = computeInsolationGrid(dt_num, latGrid, lonGrid, time_zone, rotation, dst);
            
            epsilon = 1;  % Minimal threshold to avoid division by 0
            insolation_grid( insolation_grid < epsilon) = epsilon;

            target_insol = median(insolation_grid(:));
            target_max = max(data,[],'all');

            min_insol = 15; % in W/m² 
            insol_adj = max(insolation_grid, min_insol);

            %The factor is the ratio between the target value and the local insolation.
            corr_factor = target_insol ./ insol_adj;
            % Force factor to 1 min => no “reduction” in lit areas
            corr_factor(corr_factor < 1) = 1;

            % Limit correction factor to avoid extreme corrections
            max_corr = 5; % Do not multiply by more than 5
            corr_factor(corr_factor > max_corr) = max_corr;

            % 5. Apply correction on radiance
            data = data .* corr_factor;
            data(data>target_max)=target_max;
            data = data.';    

            % (3) Process by methodName
            data_preprocessed = processDataMethod(data, methodName, ...
                clipMinHP, clipMaxHP, ...
                lpWidth20, lpWidth50, lpWidth100);

            fill_value = prctile(data_preprocessed(:), VIS_fillPercentile);
            data_preprocessed(nan_mask') = fill_value;

        otherwise
            % Fallback
            warning('Unrecognized dataType: %s. Using fallback method.', dataType);
            data_preprocessed = processDataMethod(data, methodName, ...
                clipMinHP, clipMaxHP, ...
                lpWidth20, lpWidth50, lpWidth100);
    end
end
%--------------------------------------------------------------------------

function img_processed = processDataMethod(data, methodName, ...
    clipMinHP, clipMaxHP, ...
    lpWidth20, lpWidth50, lpWidth100)
% PROCESSDATAMETHOD
%  Switch among various ways to process raw or IR/VIS data.

    switch lower(methodName)
        case 'none'
            % Just transpose for standard orientation
            img_processed = data';
            
        case 'raw_normalized'
            img_processed = normalizeData(data);
            img_processed = 1 - img_processed;
            img_processed = img_processed';

        case 'truncated'
            % Example: clamp between 280 & 292.5 K
            lower_bound = 280; 
            upper_bound = 292.5;
            img_p = data;
            img_p(img_p < lower_bound) = lower_bound;
            img_p(img_p > upper_bound) = upper_bound;
            img_p = (img_p - lower_bound) / (upper_bound - lower_bound);
            img_p = 1 - img_p;
            img_processed = img_p';

        case 'highpass_20'
            img_p = applyHighPass(data, lpWidth20, false, clipMinHP, clipMaxHP);
            img_processed = img_p';

        case 'highpass_100'
            img_p = applyHighPass(data, lpWidth100, false, clipMinHP, clipMaxHP);
            img_processed = img_p';

        case 'highpass_100_sqrt'
            img_p = applyHighPass(data, lpWidth100, true, clipMinHP, clipMaxHP);
            img_processed = img_p';

        case 'highpass_50'
            img_p = applyHighPass(data, lpWidth50, false, clipMinHP, clipMaxHP);
            img_processed = img_p';

        case 'highpass_50_sqrt'
            img_p = applyHighPass(data, lpWidth50, true, clipMinHP, clipMaxHP);
            img_processed = img_p';

        case 'highpass_20_sqrt'
            img_p = applyHighPass(data, lpWidth20, true, clipMinHP, clipMaxHP);
            img_processed = img_p';

        case 'raw_masked'
            maskLimit = 295;  % Example limit
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
%--------------------------------------------------------------------------

function img_out = applyHighPass(data, filterWidth, doSqrtEnhance, clipMinHP, clipMaxHP)
% APPLYHIGHPASS
%  Subtract a large Gaussian blur from the original to highlight small-scale
%  features. Optionally apply a sqrt contrast, then clamp to [clipMinHP, clipMaxHP].
    
    lowPass = imgaussfilt(data, filterWidth);
    highPass = data - lowPass;
    if doSqrtEnhance
        highPass = sqrt(abs(highPass)) .* sign(highPass);
    end
    highPass(highPass < clipMinHP) = clipMinHP;
    highPass(highPass > clipMaxHP) = clipMaxHP;
    img_out = (highPass - clipMinHP) / (clipMaxHP - clipMinHP);
    img_out = 1 - img_out;
end
%--------------------------------------------------------------------------

function data_win = applyRadialWindow(data_in, radius_factor, decay_rate)
    % Compute global median of the input data
    median_val = median(data_in(:));
    
    % Prepare coordinate system
    [rows, cols] = size(data_in);
    cx = cols/2;
    cy = rows/2;
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Radial distance from center
    R = sqrt((X - cx).^2 + (Y - cy).^2);
    
    % Maximum radius
    maxR = radius_factor * min(cx, cy);
    
    % Compute window values using a logistic function
    window = 1 ./ (1 + exp(decay_rate * (R - maxR)));
    
    % Instead of simply data_in .* window, blend with the median:
    %   - Where window = 1, output ~ data_in
    %   - Where window = 0, output ~ median_val
    data_win = window .* data_in + (1 - window) .* median_val;
end
%--------------------------------------------------------------------------

function data_win = applyRectangularWindow(data_in, radius_factor, decay_rate)
    % Compute global median of the input data
    median_val = median(data_in(:));
    
    % Prepare coordinate system
    [rows, cols] = size(data_in);
    cx = cols/2;
    cy = rows/2;
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Normalized absolute distances from center
    dx = abs(X - cx) / cx;
    dy = abs(Y - cy) / cy;
    
    % R is the "rectangular" distance metric, i.e. the maximum of dx, dy
    R = max(dx, dy);
    
    % Compute window values using a logistic function
    window = 1 ./ (1 + exp(decay_rate * (R - radius_factor)));
    
    % Blend data_in with the median (instead of fading to zero):
    data_win = window .* data_in + (1 - window) .* median_val;
end
%--------------------------------------------------------------------------

function produceAnnotatedImages(dataType, spec_full, data_background, squares, ...
    Scales, Angles, outDir, frameDateStr, ...
    saverose, nAngles_fineFactor, nScales_fineFactor, ...
    peakDetectionFactor,contourOption, contourArray)
% PRODUCEANNOTATEDIMAGES
%  Creates a wave-rose plot (in polar form, duplicated for ± angles),
%  detects peaks based on a threshold, and creates overlay images for each
%  peak region with wavelet real-part and power contours.


    %% 1) Dimensions & Basic Summaries
    [Ny_sh, Nx_sh, nScales, nAngles] = size(spec_full);
    [Ny_orig, Nx_orig] = size(data_background);

    power = abs(spec_full).^2;  % wavelet power
    innerpower = squeeze(mean(mean(power, 1, 'omitnan'), 2, 'omitnan'));

    %% 2) Wave-Rose Visualization with Logarithmic Radial Scale
    % 2.1) Interpolate power to a finer grid
    Angles_fine = linspace(min(Angles), max(Angles), nAngles_fineFactor*nAngles);
    % Use a logarithmic spacing for scales:
    % First, create a linear grid and then take the logarithm.
    % Alternatively, you could directly use logspace. Here we demonstrate by
    % applying the logarithm to the linear grid.
    Scales_fine_linear = linspace(min(Scales), max(Scales), nScales_fineFactor*nScales);
    
    % Create meshgrids for the original and fine grids:
    [Theta_orig, R_orig] = meshgrid(Angles, Scales);
    [Theta_fine, R_fine] = meshgrid(Angles_fine, Scales_fine_linear);
    
    % Transform the radial coordinate to a logarithmic scale.
    % (Assumes Scales > 0.)
    R_orig_log = log10(R_orig);
    R_fine_log = log10(R_fine);
    
    % Interpolate the innerpower on the original (linear) grid and then use
    % the logarithmic radial coordinates for plotting.
    F = griddedInterpolant(Theta_orig', R_orig', innerpower', 'spline');
    innerpower_fine = F(Theta_fine', R_fine')';
    
    % Compute cartesian coordinates for the fine grid using the log-transformed radius:
    [X_pos_fine, Y_pos_fine] = pol2cart(Theta_fine, R_fine_log);
    [X_neg_fine, Y_neg_fine] = pol2cart(Theta_fine + pi, R_fine_log);
    
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
    
    uistack(ax1, 'top');
    linkprop([ax1 ax2], {'XLim','YLim','Position','CameraPosition','CameraUpVector'});
    
    % 2.2) Peak detection based on threshold (for labeling)
    threshold_orig = mean(innerpower(:)) + peakDetectionFactor * std(innerpower(:));
    bwMask_orig = (innerpower >= threshold_orig);
    CC = bwconncomp(bwMask_orig, 4);
    numPeaks = CC.NumObjects;
    
    % Draw contour lines on the coarse (original) polar grid.
    [X_orig, Y_orig] = pol2cart(Theta_orig, R_orig_log);
    contour(ax1, X_orig, Y_orig, innerpower, [threshold_orig threshold_orig], 'r-', 'LineWidth',2);
    contour(ax2, X_orig, Y_orig, innerpower, [threshold_orig threshold_orig], 'r-', 'LineWidth',2);
    
    % Label each connected region (peak) on the polar plot.
    for pk = 1:numPeaks
        [scaleIndices, angleIndices] = ind2sub(size(bwMask_orig), CC.PixelIdxList{pk});
        meanScale = mean(Scales(scaleIndices), 'omitnan');
        meanAngle = mean(Angles(angleIndices), 'omitnan');
        [x_peak, y_peak] = pol2cart(meanAngle, log10(meanScale));
        text(ax1, x_peak, y_peak, sprintf('%d', pk), ...
            'Color','k','FontWeight','bold', ...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    
    %% 2.3) Annotations and Grid Overlays (using logarithmic radial axis)
    % Radial circles: plot circles at each original scale, but use log10(scales)
    for i = 1:length(Scales)
        theta_ring = linspace(0, 2*pi, 100);
        [x_ring, y_ring] = pol2cart(theta_ring, log10(Scales(i)));
        plot(ax1, x_ring, y_ring, 'k--', 'LineWidth',0.5);
        plot(ax2, x_ring, y_ring, 'k--', 'LineWidth',0.5);
        % Label the circle with the original scale value.
        text(ax1, log10(Scales(i))*1.05, 0, sprintf('%.1f', Scales(i)), ...
            'HorizontalAlignment','left','FontSize',8);
    end
    
    % Angular lines
    angle_ticks = linspace(0, 2*pi, 13);
    angle_labels = {'0','\pi/6','\pi/3','\pi/2','2\pi/3','5\pi/6','\pi',...
                    '7\pi/6','4\pi/3','3\pi/2','5\pi/3','11\pi/6','2\pi'};
    max_r = log10(max(Scales)) * 1.1;
    
    for i = 1:length(angle_ticks)
        [x_label, y_label] = pol2cart(angle_ticks(i), max_r);
        line(ax1, [0 x_label], [0 y_label], 'Color',[0.5 0.5 0.5],'LineStyle','--');
        line(ax2, [0 x_label], [0 y_label], 'Color',[0.5 0.5 0.5],'LineStyle','--');
        if angle_ticks(i) <= pi
            text(ax1, x_label*1.05, y_label*1.05, angle_labels{i}, ...
                'HorizontalAlignment','center','FontSize',8);
        else
            text(ax2, x_label*1.05, y_label*1.05, angle_labels{i}, ...
                'HorizontalAlignment','center','FontSize',8);
        end
    end
    
    % Colorbar for the polar plot
    c = colorbar(ax1, 'Location','eastoutside');
    c.Label.String = 'Wavelet Power';
    c.Label.FontWeight = 'bold';
    ax1_pos = ax1.Position;
    ax1_pos(3) = ax1_pos(3) * 0.85;
    ax1.Position = ax1_pos;
    ax2.Position = ax1_pos;
    
    title(ax1, sprintf('Polar Wave-Rose: %d Significant Regions', numPeaks), ...
        'FontSize',12, 'FontWeight','bold');
    
    if saverose
        roseName = fullfile(outDir, sprintf('WaveRose_%s.png', frameDateStr));
        exportgraphics(figRose, roseName, 'Resolution',300);
    end
    close(figRose);
    
    %% 3) Region Summaries & Overlays (Final Annotated Image)
    scaleFactorX = Nx_orig / Nx_sh;
    scaleFactorY = Ny_orig / Ny_sh;
    
    peakRegions = cell(numPeaks,1);
    for pk = 1:numPeaks
        [scaleIndices, angleIndices] = ind2sub(size(bwMask_orig), CC.PixelIdxList{pk});
        scales     = Scales(scaleIndices);
        angles_deg = rad2deg(Angles(angleIndices));
    
        scale_str = join(split(num2str(scales,'%.1f ')), '/');
        angle_str = join(split(num2str(angles_deg,'%.0f ')), '/');
    
        peakRegions{pk} = struct('ScaleIndices',scaleIndices,...
                                 'AngleIndices',angleIndices,...
                                 'ScaleStr',scale_str{1},...
                                 'AngleStr',angle_str{1});
    end
    
    for pk = 1:numPeaks
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
    
        waveSum_up   = imresize(waveSum,   [Ny_orig, Nx_orig]);
        wavePower_up = imresize(wavePower, [Ny_orig, Nx_orig]);
    
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
    
        % ----- NEW CONTOUR LEVEL SYSTEM -----
        % Choose contour levels based on either absolute values or percentiles.
        switch lower(contourOption)
            case 'absolute'
                % Use the provided absolute values (assumed positive) for contours.
                % Draw red contours at the positive levels and blue contours at the corresponding negative levels.
                contourLevels = contourArray;
            case 'percentile'
                % Compute the given percentiles on the absolute values of waveSum_up.
                contourLevels = prctile(abs(waveSum_up(:)), contourArray);
            otherwise
                error('Unknown contour option. Choose either "absolute" or "percentile".');
        end
    
        % Plot contours:
        % For positive values:
        contour(waveSum_up, contourLevels, 'LineColor','red', 'LineWidth',0.5);
        % For negative values (mirror the levels):
        contour(waveSum_up, -contourLevels, 'LineColor','blue', 'LineWidth',0.5);
    
        % ----- Draw ROI squares in final image -----
        for sq = 1:numel(squares)
            xPos_orig = squares(sq).x_range(1) * scaleFactorX;
            yPos_orig = squares(sq).y_range(1) * scaleFactorY;
            w_orig    = length(squares(sq).x_range) * scaleFactorX;
            h_orig    = length(squares(sq).y_range) * scaleFactorY;
            rectangle('Position',[xPos_orig, yPos_orig, w_orig, h_orig],...
                'EdgeColor','k','LineWidth',1);
        end
    
        titleText = {sprintf('Instrument X - %s', frameDateStr), ...
                     sprintf('Peak %d/%d - Scales: %s', pk, numPeaks, peakRegions{pk}.ScaleStr), ...
                     sprintf('Angles: %s°', peakRegions{pk}.AngleStr)};
        title(titleText, 'Color','k','FontWeight','bold','FontSize',10,'Interpreter','none');
    
        outName = fullfile(outDir, sprintf('Frame_%s_Region%02d.png', frameDateStr, pk));
        saveas(fig, outName);
        close(fig);
    end
end
%--------------------------------------------------------------------------

function renameAndOrganizeFiles(dataType, startDate, endDate, sourceRootDir, outRootDir)
% RENAMEANDORGANIZEFILES
%  Scans source directory, searching for .nc files. Extracts their timestamps,
%  checks if they are in [startDate, endDate], rounds to nearest quarter-hour,
%  renames them to "INSTRUMENT_YYYY_MM_DD_HH_MM.nc" and copies them to:
%     outRootDir\<dataType>\Data

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

    destDir = fullfile(outRootDir, upper(dataType), 'Data');
    if ~exist(destDir, 'dir')
        mkdir(destDir);
    end

    ncFiles = dir(fullfile(originalDir, '*.nc*'));
    if isempty(ncFiles)
        warning('No .nc files found in %s', originalDir);
        return;
    end

    for i = 1:numel(ncFiles)
        oldName = ncFiles(i).name;
        oldPath = fullfile(originalDir, oldName);

        tsStr = extractBetween(oldName, '_s', '_e');
        if isempty(tsStr)
            continue;
        end
        tsStr = tsStr{1};

        try
            fileTS = datetime(tsStr, 'InputFormat', 'uuuuDDDHHmmssSSS');
        catch
            continue; 
        end

        if fileTS < startDate || fileTS > endDate
            continue;
        end

        fileTS_rounded = roundToQuarterHour(fileTS);
        YYYY = year(fileTS_rounded);
        MM   = month(fileTS_rounded);
        DD   = day(fileTS_rounded);
        HH   = hour(fileTS_rounded);
        MN   = minute(fileTS_rounded);

        newName = sprintf('%s_%04d_%02d_%02d_%02d_%02d.nc', ...
                          upper(dataType), YYYY, MM, DD, HH, MN);

        newPath = fullfile(destDir, newName);
        if ~isfile(newPath)
            copyfile(oldPath, newPath);
            fprintf('Copied: %s -> %s\n', oldName, newName);
        else
            fprintf('[!] File %s already exists. Skipping.\n', newName);
        end
    end
end
%--------------------------------------------------------------------------

function roundedDT = roundToQuarterHour(originalDT)
% ROUNDTOQUARTERHOUR
%  Rounds a datetime object to the nearest 15-minute interval.

    minutesFromQuarter = mod(minute(originalDT), 15);
    if minutesFromQuarter < 7.5
        roundedDT = dateshift(originalDT, 'start', 'minute') - minutes(minutesFromQuarter);
    else
        roundedDT = dateshift(originalDT, 'start', 'minute') + minutes(15 - minutesFromQuarter);
    end
end
%--------------------------------------------------------------------------

function out = normalizeData(data)
% NORMALIZEDATA
%  Scales data to [0,1] across its range, ignoring NaNs.

    mn = min(data(:));
    mx = max(data(:));
    out = (data - mn) / (mx - mn);
end
%--------------------------------------------------------------------------

function out = normalizeDataNaN(data)
% NORMALIZEDATANAN
%  Scales data to [0,1], ignoring NaNs (filling them with min or zero).

    nanMask = isnan(data);
    data(nanMask) = min(data(~nanMask));
    out = normalizeData(data);
end
%--------------------------------------------------------------------------

function insolation = computeInsolationGrid(datetime_val, latGrid, lonGrid, time_zone, rotation, dst)
% COMPUTEINSOLATIONGRID_PARALLEL Computes solar insolation (W/m²) over a latitude/longitude grid
% using parallel computing.
%
%   insolation = computeInsolationGrid_parallel(datetime_val, latGrid, lonGrid, time_zone, rotation, dst)
%
%   Inputs:
%       datetime_val : Acquisition time in datenum format
%       latGrid      : Latitude matrix (degrees)
%       lonGrid      : Longitude matrix (degrees)
%       time_zone    : Time zone offset (hours)
%       rotation     : System rotation (degrees)
%       dst          : Daylight saving time flag (true/false)
%
%   Output:
%       insolation   : Matrix of the same size as latGrid, containing insolation in W/m²
%
% The model used is a simple approximation:
%   I = I0 * cosd(zenith)  for daytime pixels (zenith < 90°)
%   I = 0                  for nighttime pixels

    % Solar constant
    I0 = 1367; % W/m²

    [nRows, nCols] = size(latGrid);
    insolation = zeros(nRows, nCols);

    % Parallelized loop processing each row independently
    parfor i = 1:nRows
        % Temporary row storage
        tempRow = zeros(1, nCols);
        for j = 1:nCols
            % Compute solar position for pixel (i, j)
            [angles, ~] = solarPosition(datetime_val, latGrid(i,j), lonGrid(i,j), time_zone, rotation, dst);
            zenith = angles(1);  % Zenith angle in degrees

            if zenith < 90
                tempRow(j) = I0 * cosd(zenith);
            else
                tempRow(j) = 0; % Nighttime pixel
            end
        end
        % Assign the computed row back to the output matrix
        insolation(i, :) = tempRow;
    end

end
%--------------------------------------------------------------------------

function [angles,projection] = solarPosition(datetime,latitude,longitude, ...
                                             time_zone,rotation,dst)
%SOLARPOSITION Calculate solar position using most basic algorithm
%   This is the most basic algorithm. It is documented in Seinfeld &
%   Pandis, Duffie & Beckman and Wikipedia.
%
% [ANGLES,PROJECTION] = SOLARPOSITION(DATE,TIME,LATITUDE,LONGITUDE,TIME_ZONE)
% returns ZENITH & AZIMUTH for all DATE & TIME pairs at LATITUDE, LONGITUDE.
% ANGLES = [ZENITH,AZIMUTH] and PROJECTION = [PHI_X, PHI_Y]
% PHI_X is projection on x-z plane & PHI_Y is projection on y-z plane.
% DATETIME can be string, vector [YEAR, MONTH, DAY, HOURS, MINUTES, SECONDS],
%   cellstring or matrix N x [YEAR, MONTH, DAY, HOURS, MINUTES, SECONDS] for N
%   times.
% LATITUDE [degrees] and LONGITUDE [degrees] are the coordinates of the site.
% TIME_ZONE [hours] of the site.
% ROTATION [degrees] clockwise rotation of system relative to north.
% DST [logical] flag for daylight savings time, typ. from March to November
%   in the northern hemisphere.
%
% References:
% http://en.wikipedia.org/wiki/Solar_azimuth_angle
% http://en.wikipedia.org/wiki/Solar_elevation_angle
%
% Mark A. Mikofski
% Copyright (c) 2013
%
%% datetime
if iscellstr(datetime) || ~isvector(datetime)
    datetime = datenum(datetime); % [days] dates & times
else
    datetime = datetime(:); % convert datenums to row
end
date = floor(datetime); % [days]
[year,~,~] = datevec(date);
time = datetime - date; % [days]
%% constants
toRadians = @(x)x*pi/180; % convert degrees to radians
toDegrees = @(x)x*180/pi; % convert radians to degrees
%% Equation of time
d_n = mod(date-datenum(year,1,1)+1,365); % day number
B = 2*pi*(d_n-81)/365; % ET parameter
ET = 9.87*sin(2*B)-7.53*cos(B)-1.5*sin(B); % [minutes] equation of time
% approximate solar time
solarTime = ((time*24-double(dst))*60+4*(longitude-time_zone*15)+ET)/60/24;
latitude_rad = toRadians(latitude); % [radians] latitude
rotation_rad = toRadians(rotation); % [radians] field rotation
t_h = (solarTime*24-12)*15; % [degrees] hour angle
t_h_rad = toRadians(t_h); % [radians]
delta = -23.45 * cos(2*pi*(d_n+10)/365); % [degrees] declination
delta_rad = toRadians(delta); % [radians]
theta_rad = acos(sin(latitude_rad)*sin(delta_rad)+ ...
    cos(latitude_rad)*cos(delta_rad).*cos(t_h_rad)); % [radians] zenith
theta = toDegrees(theta_rad); % [degrees] zenith
elevation = 90 - theta; % elevation
day = elevation>0; % day or night?
cos_phi = (cos(theta_rad)*sin(latitude_rad)- ...
    sin(delta_rad))./(sin(theta_rad)*cos(latitude_rad)); % cosine(azimuth)
% azimuth [0, 180], absolute value measured from due south, so east = west = 90,
% south = 0, north = 180
phi_south = acos(min(1,max(-1,cos_phi)));
% azimuth [0, 360], measured clockwise from due north, so east = 90,
% south = 180, and west = 270 degrees
phi_rad = NaN(size(phi_south)); % night azimuth is NaN
% shift from ATAN to ATAN2, IE: use domain from 0 to 360 degrees instead of
% from -180 to 180
phi_rad(day) = pi + sign(t_h(day)).*phi_south(day); % Shift domain to 0-360 deg
% projection of sun angle on x-z plane, measured from z-direction (up)
phi_x = toDegrees(atan2(sin(phi_rad-rotation_rad).*sin(theta_rad), ...
    cos(theta_rad))); % [degrees]
% projection of sun angle on y-z plane, measured from z-direction (up)
phi_y = toDegrees(atan2(cos(phi_rad-rotation_rad).*sin(theta_rad), ...
    cos(theta_rad))); % [degrees]
phi = toDegrees(phi_rad); % [degrees] azimuth
angles = [theta, phi]; % [degrees] zenith, azimuth
projection = [phi_x,phi_y]; % [degrees] x-z plane, y-z plane
end

%==========================================================================