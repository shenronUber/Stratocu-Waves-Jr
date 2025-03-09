% SINGLEINSTRUMENT_CROSSWAVELET_EXPLICIT
%
% This function processes a single GOES satellite instrument (e.g., IR)
% by computing single–frame wavelet transforms and then computing cross–temporal
% cross–wavelets (i.e., the coherence between consecutive frames of the same instrument).
% Only the current and previous frame's transforms are held in memory.
%
% All parameters and thresholds (for spatial scaling, wavelet analysis,
% windowing, preprocessing, etc.) are defined in the "Variables Setup" section.

%% 1) VARIABLES SETUP
%----------- DATE/TIME SETTINGS ---------------------------
startDate = datetime(2023, 10, 12, 1, 0, 0); % start processing time
endDate   = datetime(2023, 10, 12, 3, 30, 0); % end processing time

%----------- FOLDER/PATH SETTINGS ---------------------------
rootSepacDir = 'C:\Users\admin\Box\GWaves_2023_10_11-14_SEPAC';
sourceRoot   = 'C:\Users\admin\Box\GOES2go_satellite_downloads';  % (Used in renaming)

%----------- INSTRUMENT SETTINGS ----------------------------
instrument = 'IR';  % Choose 'IR' or 'VIS'

%----------- SPATIAL SCALING & RESIZING --------------------
degrees_per_pixel = 0.04;     % Degrees per pixel (typical for GOES)
km_per_degree = 111.32;       % km per degree
shrinkfactor = 2;             % Image is resized by this factor (2 => half the resolution)
invshrinkfactor = 1 / shrinkfactor;
original_px_km = degrees_per_pixel * km_per_degree;
pixel_size_km = original_px_km * shrinkfactor;

%----------- WAVELET PARAMETERS -----------------------------
Angles = 0 : pi/(7*2) : pi;           % Wavelet angles (in radians)
Scales = [2, 4, 8, 16, 32, 64, 128];       % Wavelet scales in pixel units
NANGLES = numel(Angles);                  % Number of angles
NSCALES = numel(Scales);                  % Number of scales

CustomWavelet = false;  % This flag indicates whether to use a custom (elliptical) wavelet instead of the default built-in one.
coneAngle = pi/6;      % The variable coneAngle specifies the angular extent of the directional mask and influences the wavelet's sensitivity to orientation.
sigmaX = 0.05;         % The parameter sigmaX defines the decay rate of the wavelet's envelope along the horizontal frequency axis (ωX), affecting its horizontal resolution.
sigmaY = 1.95;         % The parameter sigmaY defines the decay rate of the wavelet's envelope along the vertical frequency axis (ωY), affecting its vertical resolution.
alpha = 0.5;           % The variable alpha is an overall radial decay factor that controls how sharply the wavelet decays in the frequency domain, thereby influencing the scale (frequency) resolution independently of the directional parameters.

%----------- WINDOWING SETTINGS -----------------------------
doWindow = true;              % Flag to apply windowing
windowType = 'rectangular';   % 'radial' or 'rectangular'
radius_factor = 0.8;          % Parameter for window function (if used)
decay_rate = 10;              % Controls steepness of window edge

%----------- SQUARE-PARTITIONING PARAMETERS ----------------
window_buffer = 0;          % Number of pixels to ignore at each edge
square_size_deg = 5;        % Square size (in degrees) for ROI partitioning

%----------- PREPROCESSING THRESHOLDS -----------------------
% For IR
IR_threshold = 277;         % IR threshold: values below are set to NaN
IR_fillPercentile = 50;     % Fill IR masked pixels with this percentile
% For VIS (not used if processing only IR, but kept for consistency)
VIS_lowerPercentile = 10;   % VIS lower bound (percentile)
VIS_upperPercentile = 99;   % VIS upper bound (percentile)
VIS_fillPercentile = 50;    % Fill VIS NaN pixels with this percentile
Insolation_Correction = false; % activate the insolation grid calculation and correction

%----------- HIGH-PASS FILTER SETTINGS ----------------------
clipMinHP = -3;             % Minimum value to clip after highpass filtering
clipMaxHP = 3;              % Maximum value to clip after highpass filtering
lowPassFilterWidth_20 = 20;
lowPassFilterWidth_50 = 50;
lowPassFilterWidth_100 = 100;

%----------- PREPROCESSING METHOD SELECTION -----------------
switch upper(instrument)
    case 'IR'
        methodName = 'highpass_50_sqrt';
    case 'VIS'
        methodName = 'none';
    otherwise
        error('Unknown instrument: %s', instrument);
end

%----------- WAVE-ROSE & PEAK DETECTION ---------------------
nAngles_fineFactor  = 4;            % Factor to refine angular resolution in the rose plot
nScales_fineFactor  = 4;            % Factor to refine scale resolution in the rose plot
peakDetectionFactor = 1;            % Threshold factor (mean + factor*std) for peak detection
contourArray        = [95 97 99];    % [Used as either percentiles or absolute values for contouring]
ArrayMode           = 'percentile'; % 'percentile' or 'absolute'

%----------- IMAGE ANNOTATIONS & OUTPUT ---------------------
saverose = true;            % Flag to save the wave–rose image
DisplayValuePower = 4*10^-3.5;
DisplayValueCoherence = 1;
DisplayValuePhase = [-pi/2;pi/2];
DisplayValueSpeed = [-30;30];

%----------- SYNTHETIC DATA SETTINGS ----------------------
syntheticWaveMode = true;   % If true, superimpose a synthetic wave on a fixed base image
driftMode         = true;   % If true, apply a drift shift each frame using circshift

% Drift parameters (in m/s) rather than pixels/frame:
drift_speed_m_s = 10;       % e.g. 10 m/s
driftAngleDeg   = 45+90;       % e.g. 45 degrees (0 = right, 90 = up) // the image is inverted !!

% Parameters for the synthetic wave:
cphase = 15;                   % Phase speed (m/s)
wavelength = 150e3;            % Wavelength in meters
direction = 235;               % Propagation direction in degrees
zamplitude = 100;              % Vertical amplitude (m)
PBLdepth = 1000;               % Boundary layer depth (m)
dB_dzPBL = 0.1;                  % dB / (dZ/PBLdepth) | Change of brightness with zamplitude


% Parameters for the spatial amplitude window (wave packet)
packet_center_x = -400e3;
packet_center_y = -400e3;
packet_width_x = 400e3;
packet_width_y = 300e3;

% Synthetic wave scaling factor for spatial coordinates
% This factor modifies DX, which represents the real-world distance per pixel.
% A smaller DX means that each pixel covers a smaller physical distance, 
% effectively "zooming in" and making the wave pattern appear larger in the image.
% Conversely, a larger DX means that each pixel covers a larger real-world distance, 
% making the wave pattern appear *smaller* and more compressed.
%
% Example: 
%   - DXFactor = 1 means the default scaling (1:1 with pixel size).
%   - DXFactor = 1/4 means the wave features are stretched, appearing 4x larger.
%   - DXFactor = 4 means the wave features are shrunk, appearing 4x smaller.
DXFactor = 1/4;  

% Seconds between frames (important for drift or wave stepping)
time_resolution = 1800;

%% 2) RETRIEVE FILE LIST
% Raw data is assumed to be in:
%    <rootSepacDir>\INSTRUMENT\Data
dataDir = fullfile(rootSepacDir, upper(instrument), 'Data');
if ~exist(dataDir, 'dir')
    error('Data directory for %s not found: %s', instrument, dataDir);
end
[fNames, fTimes, varName] = getDateRangeFiles(dataDir, startDate, endDate);
numFrames = numel(fTimes);
if numFrames == 0
    fprintf('No frames found for %s in the given time period.\n', instrument);
    return;
end
fprintf('Found %d frames for instrument %s.\n', numFrames, instrument);

%% 3) MAIN PROCESSING LOOP (Single instrument + cross–temporal coherence)
%    Rewritten to accumulate sums over the entire domain first.

% 3A) Allocate accumulators for single–frame wave–rose (i.e., power) 
%     and for cross–temporal coherence. We accumulate across all frames.
[rowsF, colsF] = deal([]);  % Will be set after first frame is processed
power_sum = [];             % For sum of |spec_full|^2  across frames
crossSpec_sum = [];         % For sum of (B1 * conj(B2)) across frame pairs
B1_auto_sum   = [];         % For sum of |B1|^2
B2_auto_sum   = [];         % For sum of |B2|^2
phaseExp_sum = [];

prevWaveletSpec = [];  % Will store wavelet from previous frame for cross–temporal coherence

for f_idx = 1:numFrames
    % ---- Prepare basic info about this frame ----
    thisTime = fTimes(f_idx);
    frameDateStr = datestr(thisTime, 'yyyy_mm_dd_HHMMSS');
    fprintf('\nProcessing frame [%d/%d]: %s\n', f_idx, numFrames, frameDateStr);

    % (Optional) define an output dir per frame, 
    % but we’re no longer doing per‐frame ROI wave–rose.
    singleOutDir = fullfile(rootSepacDir, 'Test');
    if ~exist(singleOutDir, 'dir')
        mkdir(singleOutDir);
    end
    singleNcFile = fullfile(singleOutDir, sprintf('FrameWavelet_%s.nc', frameDateStr));

    % ---- Read the raw data from file ----
    thisFileName = fNames{f_idx};
    thisFullPath = fullfile(dataDir, thisFileName);
    data = double(ncread(thisFullPath, varName));

    % ---- On the first frame, store base_frame & possibly precompute synthetic-wave grids ----
    if f_idx == 1
        base_frame = data;

        if syntheticWaveMode
            [rowsF, colsF] = size(base_frame);
            [X, Y] = meshgrid(1:colsF, 1:rowsF);

            % Real‐world pixel size (m), factoring in DXFactor:
            DX = 1000 * pixel_size_km * DXFactor;
            Xm = (X - mean(X(:))) * DX;
            Ym = (Y - mean(Y(:))) * DX * -1;  % negative if Y runs downward
        end

    elseif syntheticWaveMode || driftMode
        % If doing synthetic wave or drift, start from the same base frame each iteration:
        data = base_frame;
    end

    % ---- Synthetic wave injection (if enabled) ----
    if syntheticWaveMode
        % Compute wave parameters
        k = (2 * pi / wavelength) * cosd(direction);
        l = (2 * pi / wavelength) * sind(direction);
        omega = cphase * (2 * pi / wavelength);
        
        % Time for current frame
        t = (f_idx - 1) * time_resolution;  % e.g. in seconds
        
        % Evolving phase
        phase = k * Xm + l * Ym - omega * t;
        
        % Vertical displacement
        dz = zamplitude * sin(phase);
        
        % Envelope to localize wave in a region
        Ampwindow = exp( -(((Xm - packet_center_x) / packet_width_x).^2 ...
                        + ((Ym - packet_center_y) / packet_width_y).^2) );
        dz = dz .* Ampwindow;

        % Optionally, compute horizontal displacements:
        dxy = (zamplitude / PBLdepth) * wavelength * sin(phase - pi/2) / DX;
        dx = dxy .* cosd(direction) .* Ampwindow;
        dy = dxy .* sind(direction) .* Ampwindow;
        
        % Create new coordinates for interpolation:
        XI = X - dx;
        YI = Y - dy;

        % Warp the fixed base image:
        warped_img = interp2(X, Y, base_frame, XI, YI, 'linear', 0);
        % Modulate brightness with the vertical displacement (dz):
        modulated_img = warped_img .* (1 + dz / PBLdepth * dB_dzPBL);
        
        % Use the resulting image as the data for further processing:
        data = modulated_img;
    end

    % ---- Drift shift (if enabled) ----
    if driftMode && f_idx > 1
        % Determine how many meters the image should shift in the time between frames:
        driftDistance_m = drift_speed_m_s * time_resolution* (f_idx - 1);  % e.g. 10 m/s * 1800 s = 18000 m
        
        % Convert to pixel shift
        driftDistance_km = driftDistance_m / 1000;
        pxShift = driftDistance_km / (pixel_size_km);  % e.g. if pixel_size_km=9 => pxShift=2000/9
        pxShift = round(pxShift);  % round to nearest integer for circshift
        
        % Break it into X shift and Y shift based on the drift angle
        % (Angle=0 => shift in +X direction, 90 => shift in -Y direction, etc.)
        shift_dx =  pxShift * cosd(driftAngleDeg);
        shift_dy = -pxShift * sind(driftAngleDeg); 
        % The minus sign depends on how you define the Y orientation in your matrix:
        % Typically, "down" in the matrix is +Y, so if angle=90 means shift "up", that's negative row shift.
        
        shift_dx = round(shift_dx);
        shift_dy = round(shift_dy);
        
        % circshift( data, [row_shift, column_shift] )
        data = circshift(data, [shift_dy, shift_dx]);
    end
%%
    % ---- Preprocessing (thresholds, highpass, etc.) ----
    data_pre = preprocessFrame(data, instrument, methodName, ...
                thisFullPath, thisTime, ...
                IR_threshold, IR_fillPercentile, ...
                VIS_lowerPercentile, VIS_upperPercentile, ...
                VIS_fillPercentile, clipMinHP, clipMaxHP, ...
                lowPassFilterWidth_20, lowPassFilterWidth_50, lowPassFilterWidth_100, ...
                Insolation_Correction);
%%
    % ---- Resize or window if needed ----
    if shrinkfactor ~= 1
        data_pre = imresize(data_pre, invshrinkfactor);
    end
    if doWindow
        switch lower(windowType)
            case 'radial'
                data_pre = applyRadialWindow(data_pre, radius_factor, decay_rate);
            case 'rectangular'
                data_pre = applyRectangularWindow(data_pre, radius_factor, decay_rate);
        end
    end

    % Determine the final 2D size (post‐resize):
    [rowsF, colsF] = size(data_pre);

    % ---- Compute wavelet transform on this processed frame ----
    if CustomWavelet
        waveStruct = barebonesCauchy2D_Elliptical_NoShift( ...
                         data_pre, Scales, Angles, ...
                         coneAngle, sigmaX, sigmaY, alpha);
    else
        waveStruct = cwtft2(data_pre, ...
                            'wavelet','cauchy', ...
                            'scales', Scales, ...
                            'angles', Angles);
    end

    % cfs: size [rowsF, colsF, NSCALES, NANGLES]
    spec_full = squeeze(waveStruct.cfs);

    % Optional amplitude scaling by 2/scale:
    for iS = 1:NSCALES
        spec_full(:,:,iS,:) = spec_full(:,:,iS,:) * (2 / Scales(iS));
    end

    % ---- Init accumulators if this is the first time we know rowsF, colsF ----
    if isempty(power_sum)
        power_sum      = zeros(rowsF, colsF, NSCALES, NANGLES, 'like', spec_full);
        crossSpec_sum  = zeros(rowsF, colsF, NSCALES, NANGLES, 'like', spec_full);
        B1_auto_sum    = zeros(rowsF, colsF, NSCALES, NANGLES, 'like', spec_full);
        B2_auto_sum    = zeros(rowsF, colsF, NSCALES, NANGLES, 'like', spec_full);

        % For phase difference averaging:
        phaseExp_sum   = zeros(rowsF, colsF, NSCALES, NANGLES, 'like', spec_full) + 0i;
    end

    % ---- Accumulate single‐frame wave–rose (power) sums ----
    power_sum = power_sum + abs(spec_full).^2;

    %--- Cross-temporal stuff only if we have a previous frame ---
    if ~isempty(prevWaveletSpec)
        crossSpec_product = prevWaveletSpec .* conj(spec_full);
        crossSpec_sum = crossSpec_sum + crossSpec_product;

        B1_auto_sum = B1_auto_sum + abs(prevWaveletSpec).^2;
        B2_auto_sum = B2_auto_sum + abs(spec_full).^2;

        % Also accumulate a sum of the phase difference:
        phase_mat = angle(crossSpec_product);
        phaseExp_sum = phaseExp_sum + exp(1i * phase_mat);
    end

    prevWaveletSpec = spec_full;  % Store for next iteration

end  % end for f_idx=1:numFrames

%% 4A) ROI-BASED SUMMARIES

% Build ROI squares.

x_buffer_range = (window_buffer+1):(colsF-window_buffer);
y_buffer_range = (window_buffer+1):(rowsF-window_buffer);
adjusted_frame_width = length(x_buffer_range);
adjusted_frame_height = length(y_buffer_range);
square_size_px = round(square_size_deg / degrees_per_pixel);
num_squares_x = ceil(adjusted_frame_width / square_size_px);
num_squares_y = ceil(adjusted_frame_height / square_size_px);
squares = [];
idxS = 1;
for iy = 1:num_squares_y
    for ix = 1:num_squares_x
        x_start = floor((ix - 1) * adjusted_frame_width / num_squares_x) + 1;
        y_start = floor((iy - 1) * adjusted_frame_height / num_squares_y) + 1;
        x_end = floor(ix * adjusted_frame_width / num_squares_x);
        y_end = floor(iy * adjusted_frame_height / num_squares_y);
        if x_end > x_start && y_end > y_start
            squares(idxS).x_range = x_buffer_range(x_start:x_end);
            squares(idxS).y_range = y_buffer_range(y_start:y_end);
            squares(idxS).index = idxS;
            idxS = idxS + 1;
        end
    end
end

numSquares = numel(squares);
totalFrames = numFrames;          % For single-frame power average
numPairs    = (numFrames - 1);    % For cross-temporal pairs

% For each ROI, we’ll compute:
%   - Average single-frame power => (power_sum / totalFrames)
%   - Coherence => crossSpec_sum / autoSpec_sums
%   - Phase average => (phaseDiff_sum / numPairs)

numSquares   = numel(squares);
roiPowerCell = cell(numSquares,1);      % store wave–rose for power
roiCohCell   = cell(numSquares,1);      % store wave–rose for coherence
roiSpeedCell = cell(numSquares,1);      % store wave–rose for speed (or phase)

for iROI = 1 : numSquares
    xR = squares(iROI).x_range;  % e.g. [x_start : x_end]
    yR = squares(iROI).y_range;  % e.g. [y_start : y_end]

    %--------------------------
    % (A) Single-frame average power 
    %--------------------------
    localPow = power_sum(yR, xR, :, :);       % subarray
    sumPow   = sum(sum(localPow, 1, 'omitnan'), 2, 'omitnan'); 
    % sumPow => size [1,1,NSCALES,NANGLES], so squeeze:
    sumPow   = squeeze(sumPow);  % => [NSCALES, NANGLES]
    numPixROI = length(yR) * length(xR);
    avgPower = sumPow / (totalFrames * numPixROI);
   
    
    roiPowerCell{iROI} = avgPower;  % store

    %--------------------------
    % (B) Cross-temporal coherence 
    %--------------------------
    localCross = crossSpec_sum(yR, xR, :, :);
    localB1    = B1_auto_sum(yR, xR, :, :);
    localB2    = B2_auto_sum(yR, xR, :, :);

    sumCross = squeeze(sum(sum(localCross, 1, 'omitnan'), 2, 'omitnan'));  % => [NSCALES, NANGLES]
    sumB1    = squeeze(sum(sum(localB1,  1, 'omitnan'), 2, 'omitnan')); 
    sumB2    = squeeze(sum(sum(localB2,  1, 'omitnan'), 2, 'omitnan'));

    % Standard formula for coherence^2:
    % gamma^2 = |SumCross|^2 / ( SumB1 * SumB2 )
    gamma_sq = ( abs(sumCross).^2 ) ./ ( sumB1 .* sumB2 );

    roiCohCell{iROI} = gamma_sq;
    %--------------------------
    % (C) Naive average phase difference => "speed wave–rose" 
    %--------------------------
    localPhaseExp = phaseExp_sum(yR, xR, :, :);
    % Sum over the spatial dimensions (rows and columns):
    sumPhaseExp = squeeze(sum(sum(localPhaseExp, 1, 'omitnan'), 2, 'omitnan'));
    numPixROI = length(yR) * length(xR);
    avgPhase = angle(sumPhaseExp / (numPairs * numPixROI));

    roiSpeedCell{iROI} = avgPhase;

end

%% 4B) Produce Overlaid Figures with Inset Wave–Roses

% 1) Choose a background image. For instance, the last preprocessed frame:
bgFrameIndex = numFrames;  % last frame
thisFileName = fNames{bgFrameIndex};
thisFullPath = fullfile(dataDir, thisFileName);
data_bg = double(ncread(thisFullPath, varName));
data_bg_pre = preprocessFrame(data, instrument, methodName, ...
                thisFullPath, thisTime, ...
                IR_threshold, IR_fillPercentile, ...
                VIS_lowerPercentile, VIS_upperPercentile, ...
                VIS_fillPercentile, clipMinHP, clipMaxHP, ...
                lowPassFilterWidth_20, lowPassFilterWidth_50, lowPassFilterWidth_100, ...
                Insolation_Correction);  % do same steps as you do in the loop

if doWindow
    switch lower(windowType)
        case 'radial'
            data_bg_pre = applyRadialWindow(data_bg_pre, radius_factor, decay_rate);
        case 'rectangular'
            data_bg_pre = applyRectangularWindow(data_bg_pre, radius_factor, decay_rate);
    end
end

produceOverlayWaveRose('Power', roiPowerCell, squares, num_squares_x, num_squares_y, ...
                        Scales, Angles, DisplayValuePower, data_bg_pre, singleOutDir, 'Global_PowerWaveRose_Overlay.png');

produceOverlayWaveRose('Coherence', roiCohCell, squares, num_squares_x, num_squares_y, ...
                        Scales, Angles, DisplayValueCoherence, data_bg_pre, singleOutDir, 'Global_CoherenceWaveRose_Overlay.png');

produceOverlayWaveRose('Phase', roiSpeedCell, squares, num_squares_x, num_squares_y, ...
                        Scales, Angles, DisplayValuePhase, data_bg_pre, singleOutDir, 'Global_PhaseWaveRose_Overlay.png');

produceOverlayWaveRose('Speed', roiSpeedCell, squares, num_squares_x, num_squares_y, ...
                        Scales, Angles, DisplayValueSpeed, data_bg_pre, singleOutDir, 'Global_SpeedWaveRose_Overlay.png');


%% 5) GLOBAL AVERAGE WAVE ROSES
% Aafter Section 3, we have the following accumulators:
%   power_sum     : sum over frames of |spec_full|^2, size [rowsF, colsF, NSCALES, NANGLES]
%   crossSpec_sum : sum over consecutive-frame pairs of (B1 .* conj(B2))
%   B1_auto_sum   : sum over pairs of |B1|^2 (previous frame)
%   B2_auto_sum   : sum over pairs of |B2|^2 (current frame)
%   phaseExp_sum  : sum over pairs of exp(1i*phase_diff), for circular averaging

% Define total numbers:
totalFrames = numFrames;         % for single-frame (power) sums
numPairs    = numFrames - 1;      % for cross-temporal quantities
numPixels   = rowsF * colsF;       % total pixels per frame

% ----- (A) Global Power Wave–Rose -----
% Average power over frames and spatial domain:
globalPower = squeeze( sum(sum(power_sum, 1, 'omitnan'), 2, 'omitnan') ) ...
              / (totalFrames * numPixels);
% Produce the global wave–rose plot:
produceAggregatedWaveRose('Power_Global', globalPower, Scales, Angles, ...
    singleOutDir, 'Global_PowerWaveRose', saverose, DisplayValuePower);

% ----- (B) Global Coherence Wave–Rose -----
globalCross = squeeze( sum(sum(crossSpec_sum, 1, 'omitnan'), 2, 'omitnan') );
globalB1    = squeeze( sum(sum(B1_auto_sum,   1, 'omitnan'), 2, 'omitnan') );
globalB2    = squeeze( sum(sum(B2_auto_sum,   1, 'omitnan'), 2, 'omitnan') );
% Compute coherence (using standard formula: |S12|^2 / (S11*S22)):
globalCoherence = ( abs(globalCross).^2 ) ./ (globalB1 .* globalB2 );
produceAggregatedWaveRose('Coherence_Global', globalCoherence, Scales, Angles, ...
    singleOutDir, 'Global_CoherenceWaveRose', saverose, DisplayValueCoherence);

% ----- (C) Global Speed (Phase) Wave–Rose -----
% For phase differences, we average the complex exponentials (for circular averaging)
globalPhaseExp = squeeze( sum(sum(phaseExp_sum, 1, 'omitnan'), 2, 'omitnan') );
% Average over the number of pairs and spatial domain:
globalAvgPhase = angle( globalPhaseExp / (numPairs * numPixels) );
produceAggregatedWaveRose('Speed_Global', globalAvgPhase, Scales, Angles, ...
    singleOutDir, 'Global_SpeedWaveRose', saverose, DisplayValueSpeed);

% ----- (D) Global Phase Wave–Rose -----
% For phase differences, we average the complex exponentials (for circular averaging)
produceAggregatedWaveRose('Phase_Global', globalAvgPhase, Scales, Angles, ...
    singleOutDir, 'Global_PhaseWaveRose', saverose, DisplayValuePhase);


fprintf('\nAll done. Single–frame and cross–temporal wavelet (coherence) computations complete.\n');


%% HELPER FUNCTIONS
%==========================================================================
function produceOverlayWaveRose(metricLabel, roiCell, squares, num_squares_x, num_squares_y, ...
                                Scales, Angles, DisplayValue, data_bg_pre, outDir, fileName)
% produceOverlayWaveRose produces a figure with insets of wave–rose plots
% overlaid on a background image.
%
% Inputs:
%   metricLabel   - A string indicating the metric (e.g. 'Power', 'Coherence', 'Speed').
%   roiCell       - A cell array containing the [nScales x nAngles] wave–rose matrix for each ROI.
%   squares       - A structure array defining each ROI with fields 'x_range' and 'y_range'.
%   num_squares_x - Number of ROI squares in the horizontal direction.
%   num_squares_y - Number of ROI squares in the vertical direction.
%   Scales, Angles- Vectors defining your wavelet scales and angles.
%   DisplayValue  - A scalar (or two-element vector for speed) to set the color axis.
%   data_bg_pre   - The background image (preprocessed) to display.
%   outDir        - Output directory where the figure will be saved.
%   fileName      - Name of the output file (e.g. 'Global_PowerWaveRose_Overlay.png').
%
% Example:
%   produceOverlayWaveRose('Power', roiPowerCell, squares, numSquares_x, numSquares_y, ...
%                           Scales, Angles, DisplayValuePower, data_bg_pre, outDir, 'Global_PowerOverlay.png');

    % Create the figure with the background image.
    fig = figure('Units','normalized','Position',[0.1 0.1 0.6 0.6]);
    axMain = axes(fig);
    imagesc(axMain, data_bg_pre);
    colormap(axMain, 'gray');
    axis(axMain, 'image');
    axis(axMain, 'off');
    title(axMain, ['Background with ' metricLabel ' Wave–Roses']);
    
    % Get main axes position and data limits.
    posMain = get(axMain, 'Position');   % [x0, y0, w, h] in figure normalized coordinates.
    xLim    = axMain.XLim;              % e.g. [1, Nx]
    yLim    = axMain.YLim;              % e.g. [1, Ny]
    
    % Set inset scaling and offsets.
    sizeFactor = 2;    % Factor to enlarge each inset.
    xOffset   = 0.035; % Horizontal offset between insets.
    yOffset   = 0.08;  % Vertical offset between insets.
    
    % Loop over all ROIs.
    for iy = 1:num_squares_y
        for ix = 1:num_squares_x
            idxS = (iy - 1) * num_squares_x + ix;  % index into squares array.
            roi = squares(idxS);
            
            % Get pixel coordinate bounds for the ROI.
            xMin = min(roi.x_range);
            xMax = max(roi.x_range);
            yMin = min(roi.y_range);
            yMax = max(roi.y_range);
            width  = xMax - xMin + 1;
            height = yMax - yMin + 1;
            
            % Convert ROI data coordinates to figure–normalized coordinates.
            xNorm = posMain(1) + ((xMin - xLim(1)) / (xLim(2) - xLim(1))) * posMain(3);
            wNorm = (width / (xLim(2) - xLim(1))) * posMain(3) * sizeFactor;
            yNorm = posMain(2) + ((yLim(2) - yMax) / (yLim(2) - yLim(1))) * posMain(4);
            hNorm = (height / (yLim(2) - yLim(1))) * posMain(4) * sizeFactor;
            
            % Adjust for spacing between insets.
            xNorm = xNorm + (ix - 1) * xOffset;
            yNorm = yNorm - (iy - 1) * yOffset;
            
            % Create an inset axes at the computed position.
            axInset = axes('Position', [xNorm+0.07, yNorm-0.08, wNorm, hNorm]);
            
            % Retrieve the corresponding wave–rose matrix.
            roiMatrix = roiCell{idxS};
            
            % Plot the wave–rose for this ROI. (DisplayAggregatedWaveRose is your helper.)
            DisplayAggregatedWaveRose([metricLabel ' Square'], roiMatrix, Scales, Angles, axInset, DisplayValue);
            
            axis(axInset, 'off');
        end
    end
    
    % Save the overlay figure.
    overlayFileName = fullfile(outDir, fileName);
    exportgraphics(fig, overlayFileName, 'Resolution', 400);
    fprintf('Saved %s wave–rose overlay to: %s\n', metricLabel, overlayFileName);
end

function DisplayAggregatedWaveRose(labelStr, waveRoseMat, Scales, Angles, axHandle, maxVal)
% DISPLAYAGGREGATEDWAVEROSE Displays an aggregated wave-rose plot in a provided axes.
%
%   DISPLAYAGGREGATEDWAVEROSE(labelStr, waveRoseMat, Scales, Angles, axHandle, maxVal)
%
%   Inputs:
%     labelStr    - A string label (e.g., 'Coherence_Square') that will be used 
%                   to set the title and determine the colormap.
%     waveRoseMat - A 2D matrix of size [nScales x nAngles] representing the aggregated
%                   power (or coherence) spectrum.
%     Scales      - A vector of scales.
%     Angles      - A vector of angles (in radians).
%     axHandle    - A handle to an existing axes in which to display the wave-rose.
%     maxVal      - (Optional) A scalar maximum value to set the color axis to [0 maxVal].
%
%   Example:
%     % Suppose axInset is an inset axes handle and maxVal is computed globally.
%     DisplayAggregatedWaveRose('Coherence_Square', coherenceMatrix, Scales, Angles, axInset, maxVal);

    % Check that axHandle is provided.
    if nargin < 5 || isempty(axHandle)
        error('An axes handle (axHandle) must be provided.');
    end
    
    if startsWith(labelStr, 'Speed', 'IgnoreCase', true) || startsWith(labelStr, 'Phase', 'IgnoreCase', true)
        minVal = maxVal(1);
        maxVal = maxVal(2);
    end

    % Set the current axes to the provided handle and hold it.
    axes(axHandle);
    hold(axHandle, 'on');

    % Determine dimensions.
    nScales = length(Scales);
    nAngles = length(Angles);
    
    % Use the provided matrix as the "innerpower" for display.
    innerpower = waveRoseMat;
    
    % Refine grid for smoother plotting.
    nAngles_fineFactor = 4;
    nScales_fineFactor = 4;
    Angles_fine = linspace(min(Angles), max(Angles), nAngles_fineFactor * nAngles);
    Scales_fine_linear = linspace(min(Scales), max(Scales), nScales_fineFactor * nScales);
    
    [Theta_orig, R_orig] = meshgrid(Angles, Scales);
    [Theta_fine, R_fine] = meshgrid(Angles_fine, Scales_fine_linear);
    R_orig_log = log10(R_orig);
    R_fine_log = log10(R_fine);
    
    % Interpolate the innerpower onto the fine grid.
    F = griddedInterpolant(Theta_orig', R_orig', innerpower', 'spline');
    innerpower_fine = F(Theta_fine', R_fine')';
  
    if startsWith(labelStr, 'Speed', 'IgnoreCase', true)
        pixel_size_km = evalin('base','pixel_size_km');
        % For each row (corresponding to a fine-scale value), multiply by that scale.
        for iRow = 1:size(innerpower_fine,1)
            innerpower_fine(iRow, :) = (pixel_size_km*1000*(innerpower_fine(iRow, :)/(2*pi)) .* Scales_fine_linear(iRow)* pi/sqrt(2))/ 1800;
        end
    end
    
    % Convert polar coordinates to Cartesian (using a log-scale for R).
    [X_pos_fine, Y_pos_fine] = pol2cart(Theta_fine, R_fine_log);
    [X_neg_fine, Y_neg_fine] = pol2cart(Theta_fine + pi, R_fine_log);
    
    % Plot the positive half of the rose.
    pcolor(axHandle, X_pos_fine, Y_pos_fine, innerpower_fine);
    shading(axHandle, 'interp');
    
    % Choose colormap based on label.
    if startsWith(labelStr, 'Speed', 'IgnoreCase', true)
        colormap(axHandle, 'jet');
    elseif startsWith(labelStr, 'Coherence', 'IgnoreCase', true)
        colormap(axHandle, 'turbo');
    elseif startsWith(labelStr, 'Power', 'IgnoreCase', true)
        colormap(axHandle, 'parula');
    elseif startsWith(labelStr, 'Phase', 'IgnoreCase', true)
        colormap(axHandle, 'hsv');
    end
    axis(axHandle, 'equal', 'tight', 'off');
    
    % Optionally set color axis limits if maxVal is provided.
    if startsWith(labelStr, 'Speed', 'IgnoreCase', true) || startsWith(labelStr, 'Phase', 'IgnoreCase', true)
            clim(axHandle, [minVal, maxVal]);
    elseif nargin >= 6 && ~isempty(maxVal)
            clim(axHandle, [0, maxVal]);
    end
    
    % Optionally add radial grid lines.
    max_scale_log = log10(max(Scales));
    for i = 1:length(Scales)
        level_log = log10(Scales(i));
        theta_ring = linspace(0, 2*pi, 180);
        [x_ring, y_ring] = pol2cart(theta_ring, level_log);
        plot(axHandle, x_ring, y_ring, 'k--', 'LineWidth',0.5);
        %text(axHandle, level_log, 0, sprintf('%.1f', Scales(i)), ...
        %     'Color','k','FontSize',4,'HorizontalAlignment','left');
    end

    
    % Add angular lines.
    angle_ticks = linspace(0, 2*pi, 7);
    max_r = max_scale_log * 1.1;
    for i = 1:length(angle_ticks)
        [x_label, y_label] = pol2cart(angle_ticks(i), max_r);
        line(axHandle, [0 x_label], [0 y_label], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
    end
end

function produceAggregatedWaveRose(labelStr, waveRoseMat, Scales, Angles, outDir, outName, saverose, maxVal)
    if nargin < 7
        saverose = true;
    end
    if nargin < 8
        maxVal = [];  % if not provided, use auto-scale
    end

    if size(maxVal,1)==2
        minVal = maxVal(1);
        maxVal = maxVal(2);
    end

    % Dimensions
    nScales = length(Scales);
    nAngles = length(Angles);
    
    % Use waveRoseMat as the "innerpower" (power spectrum display)
    innerpower = waveRoseMat;
    
    % Refine grid for smoother plotting.
    nAngles_fineFactor  = 4;
    nScales_fineFactor  = 4;
    Angles_fine = linspace(min(Angles), max(Angles), nAngles_fineFactor * nAngles);
    Scales_fine_linear = linspace(min(Scales), max(Scales), nScales_fineFactor * nScales);
    [Theta_orig, R_orig] = meshgrid(Angles, Scales);
    [Theta_fine, R_fine] = meshgrid(Angles_fine, Scales_fine_linear);
    R_orig_log = log10(R_orig);
    R_fine_log = log10(R_fine);
    
    % Interpolate the innerpower onto the fine grid.
    F = griddedInterpolant(Theta_orig', R_orig', innerpower', 'spline');
    innerpower_fine = F(Theta_fine', R_fine')';

    % --- Modification for Speed: multiply each row by its corresponding scale ---
    if startsWith(labelStr, 'Speed', 'IgnoreCase', true)
        pixel_size_km = evalin('base','pixel_size_km');
        % For each row (corresponding to a fine-scale value), multiply by that scale.
        for iRow = 1:size(innerpower_fine,1)
            innerpower_fine(iRow, :) = (pixel_size_km*1000*(innerpower_fine(iRow, :)/(2*pi)) .* Scales_fine_linear(iRow)* pi/sqrt(2))/ 1800;
        end
    end
    
    % Convert to Cartesian coordinates for plotting.
    [X_pos_fine, Y_pos_fine] = pol2cart(Theta_fine, R_fine_log);
    [X_neg_fine, Y_neg_fine] = pol2cart(Theta_fine + pi, R_fine_log);
    
    % Create figure.
    figRose = figure('visible','off');
    set(figRose, 'Position', [100 100 600 600]);
    ax1 = axes('Position',[0.1 0.1 0.75 0.75]);
    hold(ax1, 'on');
    pcolor(ax1, X_pos_fine, Y_pos_fine, innerpower_fine);
    shading(ax1, 'interp');
    if startsWith(labelStr, 'Coherence', 'IgnoreCase', true)
        colormap(ax1, 'turbo');
    elseif startsWith(labelStr, 'Speed', 'IgnoreCase', true)
        colormap(ax1,'jet')
    elseif startsWith(labelStr, 'Power', 'IgnoreCase', true)
        colormap(ax1, 'parula');
    elseif startsWith(labelStr, 'Phase', 'IgnoreCase', true)
        colormap(ax1, 'hsv');
    end

    if startsWith(labelStr, 'Speed', 'IgnoreCase', true)|| startsWith(labelStr, 'Phase', 'IgnoreCase', true)
        clim(ax1, [minVal, maxVal]);
    elseif ~isempty(maxVal)  
        clim(ax1, [0, maxVal]);
    else
        clim(ax1, [0, max(innerpower_fine(:))]);
    end
    axis(ax1, 'equal', 'tight', 'off');
    
    % Duplicate plot for negative angles.
    ax2 = axes('Position', ax1.Position, 'Color','none', 'HitTest','off');
    hold(ax2, 'on');
    pcolor(ax2, X_neg_fine, Y_neg_fine, innerpower_fine);
    shading(ax2, 'interp');
    if ~isempty(maxVal)
        clim(ax2, [0, maxVal]);
    else
        clim(ax2, [0, max(innerpower_fine(:))]);
    end
    axis(ax2, 'equal', 'tight', 'off');
    
    uistack(ax1, 'top');
    linkprop([ax1 ax2], {'XLim','YLim','Position','CameraPosition','CameraUpVector'});
    
    % Add radial grid lines (using logarithmic scale for the scales).
    max_scale_log = log10(max(Scales));
    for i = 1:length(Scales)
        level_log = log10(Scales(i));
        theta_ring = linspace(0, 2*pi, 180);
        [x_ring, y_ring] = pol2cart(theta_ring, level_log);
        plot(ax1, x_ring, y_ring, 'k--', 'LineWidth',0.5);
        plot(ax2, x_ring, y_ring, 'k--', 'LineWidth',0.5);
        val_linear = 10^(level_log);
        text(ax1, level_log, 0, sprintf('%.2f', val_linear), 'Color','k','FontSize',4,'HorizontalAlignment','left');
    end
    
    % Add angular lines.
    angle_ticks = linspace(0, 2*pi, 7);
    max_r = max_scale_log * 1.1;
    for i = 1:length(angle_ticks)
        [x_label, y_label] = pol2cart(angle_ticks(i), max_r);
        line(ax1, [0 x_label], [0 y_label], 'Color',[0.5 0.5 0.5],'LineStyle','--');
        line(ax2, [0 x_label], [0 y_label], 'Color',[0.5 0.5 0.5],'LineStyle','--');
    end
    
    title(ax1, sprintf('%s Wave–Rose', labelStr), 'FontSize',12, 'FontWeight','bold');
    
    c = colorbar(ax1, 'Location','eastoutside');
    c.Label.String = 'Wavelet Power';
    c.Label.FontWeight = 'bold';
    
    ax1_pos = ax1.Position;
    ax1_pos(3) = ax1_pos(3) * 0.85;
    ax1.Position = ax1_pos;
    ax2.Position = ax1_pos;
    
    if saverose
        roseFileName = fullfile(outDir, sprintf('%s.png', outName));
        fprintf('Saving wave–rose plot to: %s\n', roseFileName);
        exportgraphics(figRose, roseFileName, 'Resolution',300);
    end
    close(figRose);
end
%--------------------------------------------------------------------------

function results = extractWaveletFeatures(dataType, spec_full, data_background, squares, ...
    Scales, Angles, frameDateStr, ...
    nAngles_fineFactor, nScales_fineFactor, peakDetectionFactor, contourOption, contourArray)
% EXTRACTWAVELETFEATURES
%  Extracts wavelet features such as the number of peaks, scale/angle pairs,
%  and contour positions for each frame.

    %% 1) Dimensions & Basic Summaries
    [Ny_sh, Nx_sh, nScales, nAngles] = size(spec_full);
    [Ny_orig, Nx_orig] = size(data_background);

    power = abs(spec_full).^2;  % wavelet power
    innerpower = squeeze(mean(mean(power, 1, 'omitnan'), 2, 'omitnan'));

    %% 2) Wave-Rose Visualization with Logarithmic Radial Scale
    % 2.1) Interpolate power to a finer grid
    Angles_fine = linspace(min(Angles), max(Angles), nAngles_fineFactor*nAngles);
    Scales_fine_linear = linspace(min(Scales), max(Scales), nScales_fineFactor*nScales);

    % Create meshgrids for the original and fine grids:
    [Theta_orig, R_orig] = meshgrid(Angles, Scales);
    [Theta_fine, R_fine] = meshgrid(Angles_fine, Scales_fine_linear);

    % Transform the radial coordinate to a logarithmic scale.
    R_orig_log = log10(R_orig);
    R_fine_log = log10(R_fine);

    % Interpolate the innerpower on the original (linear) grid and then use
    % the logarithmic radial coordinates for plotting.
    F = griddedInterpolant(Theta_orig', R_orig', innerpower', 'spline');
    innerpower_fine = F(Theta_fine', R_fine')';

    % Compute cartesian coordinates for the fine grid using the log-transformed radius:
    [X_pos_fine, Y_pos_fine] = pol2cart(Theta_fine, R_fine_log);
    [X_neg_fine, Y_neg_fine] = pol2cart(Theta_fine + pi, R_fine_log);

    % 2.2) Peak detection based on threshold (for labeling)
    threshold_orig = mean(innerpower(:)) + peakDetectionFactor * std(innerpower(:));
    bwMask_orig = (innerpower >= threshold_orig);
    CC = bwconncomp(bwMask_orig, 4);
    numPeaks = CC.NumObjects;

    % Draw contour lines on the coarse (original) polar grid.
    [X_orig, Y_orig] = pol2cart(Theta_orig, R_orig_log);

    % Label each connected region (peak) on the polar plot.
    peakRegions = struct();
    for pk = 1:numPeaks
        [scaleIndices, angleIndices] = ind2sub(size(bwMask_orig), CC.PixelIdxList{pk});
        meanScale = mean(Scales(scaleIndices), 'omitnan');
        meanAngle = mean(Angles(angleIndices), 'omitnan');
        peakRegions(pk).ScaleIndices = scaleIndices;
        peakRegions(pk).AngleIndices = angleIndices;
        peakRegions(pk).MeanScale = meanScale;
        peakRegions(pk).MeanAngle = meanAngle;
    end

    %% 3) Region Summaries & Overlays (Final Annotated Image)
    scaleFactorX = Nx_orig / Nx_sh;
    scaleFactorY = Ny_orig / Ny_sh;

    for pk = 1:numPeaks
        waveSum = zeros(Ny_sh, Nx_sh);
        wavePower = zeros(Ny_sh, Nx_sh);
        currentRegion = peakRegions(pk);

        for jj = 1:numel(currentRegion.ScaleIndices)
            s_idx = currentRegion.ScaleIndices(jj);
            a_idx = currentRegion.AngleIndices(jj);
            coeff = spec_full(:,:,s_idx,a_idx);
            waveSum = waveSum + real(coeff);
            wavePower = wavePower + abs(coeff).^2;
        end

        waveSum_up = imresize(waveSum, [Ny_orig, Nx_orig]);
        wavePower_up = imresize(wavePower, [Ny_orig, Nx_orig]);

        % ----- NEW CONTOUR LEVEL SYSTEM -----
        % Choose contour levels based on either absolute values or percentiles.
        switch lower(contourOption)
            case 'absolute'
                % Use the provided absolute values (assumed positive) for contours.
                contourLevels = contourArray;
            case 'percentile'
                % Compute the given percentiles on the absolute values of waveSum_up.
                contourLevels = prctile(abs(waveSum_up(:)), contourArray);
            otherwise
                error('Unknown contour option. Choose either "absolute" or "percentile".');
        end

        % Store contour positions
        contourPositions = struct();
        contourPositions.Positive = contourc(waveSum_up, contourLevels);
        contourPositions.Negative = contourc(waveSum_up, -contourLevels);

        % Store results for the current peak
        peakRegions(pk).ContourPositions = contourPositions;
    end

    % Store results for the current frame
    
    results.NumPeaks = numPeaks;
    results.PeakRegions = peakRegions;
end
%--------------------------------------------------------------------------

function velocityField = calculatePIV(image1, image2, shrinkfactor)
    % Calculate PIV between two images
    [xtable, ytable, utable, vtable, ~, ~, ~] = piv_FFTmulti(...
    image1, image2, ...          % Input images
    64, ...                      % Initial interrogation area (Pass 1)
    32, ...                      % Initial step (Pass 1)
    2, ...                       % Subpixel method (2D Gaussian)
    [], ...                      % mask_inpt (no mask)
    [], ...                      % roi_inpt (full image)
    3, ...                       % passes (3 iterations)
    32, ...                      % int2 (Pass 2 IA)
    16, ...                      % int3 (Pass 3 IA)
    0, ...                       % int4 (unused)
    '*linear', ...               % imdeform method
    0, ...                       % repeat (no repeated correlation)
    1, ...                       % mask_auto (disable autocorrelation)
    1, ...                       % do_linear_correlation (enable)
    0, ...                       % do_correlation_matrices (disable)
    0, ...                       % repeat_last_pass (no)
    0.025 ...                    % delta_diff_min (convergence threshold)
);

    % Downsample the velocity vectors by the shrinkfactor
    downsample_factor = shrinkfactor;
    xtable_ds = xtable(1:downsample_factor:end, 1:downsample_factor:end);
    ytable_ds = ytable(1:downsample_factor:end, 1:downsample_factor:end);
    utable_ds = utable(1:downsample_factor:end, 1:downsample_factor:end);
    vtable_ds = vtable(1:downsample_factor:end, 1:downsample_factor:end);

    % Combine the downsampled velocity vectors into a single matrix
    velocityField = cat(3, xtable_ds, ytable_ds, utable_ds, vtable_ds);
end
%--------------------------------------------------------------------------

function [fileNames, fileTimestamps, variableName] = getDateRangeFiles(dataDir, startDate, endDate)
    ncFiles = dir(fullfile(dataDir, '*.nc'));
    if isempty(ncFiles)
        fileNames = {};
        fileTimestamps = [];
        variableName = '';
        warning('No .nc files in %s', dataDir);
        return;
    end

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
        mm = str2double(parts{3});
        dd = str2double(parts{4});
        HH = str2double(parts{5});
        minPart = erase(parts{6}, '.nc');
        MN = str2double(minPart);

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
    lpWidth20, lpWidth50, lpWidth100,Insolation_Correction)
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

            if Insolation_Correction
                
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
            end
            
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

function img_processed = processDataMethod(data, methodName, clipMinHP, clipMaxHP, lpWidth20, lpWidth50, lpWidth100)
    switch lower(methodName)
        case 'none'
            img_processed = data';
        case 'truncated'
            lower_bound = 280;
            upper_bound = 292.5;
            img_p = data;
            img_p(img_p < lower_bound) = lower_bound;
            img_p(img_p > upper_bound) = upper_bound;
            img_p = standardizeData(img_p);
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
            maskLimit = 295;
            data(data > maskLimit) = NaN;
            img_p = standardizeDataNaN(data);
            img_processed = img_p';
        case 'raw'
            img_p = standardizeDataNaN(data);
            img_processed = img_p';
        otherwise
            error('Unknown method: %s', methodName);
    end
end
%--------------------------------------------------------------------------

function img_out = applyHighPass(data, filterWidth, doSqrtEnhance, clipMinHP, clipMaxHP)
    lowPass = imgaussfilt(data, filterWidth);
    highPass = data - lowPass;
    if doSqrtEnhance
        highPass = sqrt(abs(highPass)) .* sign(highPass);
    end
    highPass(highPass < clipMinHP) = clipMinHP;
    highPass(highPass > clipMaxHP) = clipMaxHP;
    img_out = standardizeData(highPass);
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

function out = standardizeData(data)
    meanVal = mean(data(:), 'omitnan');
    stdVal  = std(data(:),  'omitnan');
    out    = (data - meanVal) ./ stdVal;
end
%--------------------------------------------------------------------------

function out = standardizeDataNaN(data)
    nanMask = isnan(data);
    data(nanMask) = min(data(~nanMask));
    out = standardizeData(data);
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
            'HorizontalAlignment','left','FontSize',4);
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
                imagesc(data_background, [-2 2])
                colormap(flipud(gray))
            case 'VIS'
                image(data_background);
                colormap(gray)
            otherwise
                error('Unknown dataType.');
        end
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
        MM = month(fileTS_rounded);
        DD = day(fileTS_rounded);
        HH = hour(fileTS_rounded);
        MN = minute(fileTS_rounded);
        newName = sprintf('%s_%04d_%02d_%02d_%02d_%02d.nc', upper(dataType), YYYY, MM, DD, HH, MN);
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
    minutesFromQuarter = mod(minute(originalDT), 15);
    if minutesFromQuarter < 7.5
        roundedDT = dateshift(originalDT, 'start', 'minute') - minutes(minutesFromQuarter);
    else
        roundedDT = dateshift(originalDT, 'start', 'minute') + minutes(15 - minutesFromQuarter);
    end
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
