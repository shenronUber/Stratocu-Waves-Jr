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
%startDate = datetime(2023, 10, 12, 1, 0, 0); % start processing time
%endDate   = datetime(2023, 10, 12, 3, 30, 0); % end processing time
startDate = datetime(2023, 10, 11, 14, 0, 0); % start processing time
endDate = datetime(2023, 10, 11, 20, 30, 0); % end processing time


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

%----------- WAVELET PARAMETERS -----------------------------
Angles = 0 : pi/(7*2) : pi;           % Wavelet angles (in radians)
Scales = [2, 4, 8, 16, 32, 64];       % Wavelet scales in pixel units
Scales_orig = Scales;                 % That will anchor the Scales values in case there's a ShrinkFactor >1
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
radius_factor = 0.6;          % Parameter for window function (if used)
decay_rate = 10;              % Controls steepness of window edge

%----------- SQUARE-PARTITIONING PARAMETERS ----------------
window_buffer = 0;          % Number of pixels to ignore at each edge
square_size_deg = 10;        % Square size (in degrees) for ROI partitioning

%----------- PREPROCESSING THRESHOLDS -----------------------
% For IR
IR_threshold = 277;         % IR threshold: values below are set to NaN
IR_fillPercentile = 50;     % Fill IR masked pixels with this percentile
% For VIS (not used if processing only IR, but kept for consistency)
VIS_lowerPercentile = 10;   % VIS lower bound (percentile)
VIS_upperPercentile = 99;   % VIS upper bound (percentile)
VIS_fillPercentile = 50;    % Fill VIS NaN pixels with this percentile
Insolation_Correction = true; % activate the insolation grid calculation and correction

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
%contourArray        = [95 97 99];    % [Used as either percentiles or absolute values for contouring]
%ArrayMode           = 'percentile'; % 'percentile' or 'absolute'

%----------- IMAGE ANNOTATIONS & OUTPUT ---------------------
saverose = true;            % Flag to save the wave–rose image
DisplayValuePower = 4*10^-2;
DisplayValueCoherence = 1;
DisplayValuePhase = [-pi;pi];
DisplayValueSpeed = [-15;15];

%----------- ADVECTION CORRECTION SETTINGS -----------------
doAdvectionEstimation = true;     % Set to true to estimate mean advection
scalesForAdvection = [4, 8, 16, 32];    % Use these pixel scales for advection estimation (central scales)
cohThreshold = 0.2;               % Only use (scale,angle) bins where coherence (or its square root) is >= 0.2
amplitudeThreshold = 0;           % (Optional) ignore bins with amplitude below this threshold

%----------- SPEED CORRECTION SETTINGS -------------------
% you can either play with the parameters or use a specific calibration Matrix
beta = 0.8;         % Weight for smoothing toward baseline trend
decayFactor = 0.55;  % Amplitude decay applied to baseline smoothing
decaySharpness = 1.2;  % Controls how sharply the baseline trend decays with scale
upperCutoff = 16;  % Upper trusted scale
%nyquistScales= Scales_orig(1:2); % Scales that are likely to be hit by nyquist issue
nyquistScales= Scales_orig(:,Scales_orig<=8); % Scales that are likely to be hit by nyquist issue

matrixMode = true;
%calibrationMatrix = [0.5;0.5;0.5;0.4;0.335;0.17;0.083]; % Empirical calibration values
calibrationMatrix = [
     2,   0.50;
     4,   0.50;
     8,   0.50;
    16,   0.40;
    32,   0.335;
    64,   0.17;
   128,   0.083;
];

%----------- SYNTHETIC DATA SETTINGS ----------------------
syntheticWaveMode = true;   % If true, superimpose a synthetic wave on a fixed base image
driftMode         = true;   % If true, apply a drift shift each frame using circshift

% Drift parameters (in m/s) rather than pixels/frame:
drift_speed_m_s = 15;       % e.g. 10 m/s
driftAngleDeg   = 45+90;       % e.g. 45 degrees (0 = right, 90 = up) // the image is inverted !!


% Parameters for the synthetic wave:
cphase = 15;                   % Phase speed (m/s)
wavelength = 150e3;            % Wavelength in meters
direction = 235;               % Propagation from direction in degrees
zamplitude = 100;              % Vertical amplitude (m)
PBLdepth = 1000;               % Boundary layer depth (m)
dB_dzPBL = 0.1;                  % dB / (dZ/PBLdepth) | Change of brightness with zamplitude


% Parameters for the spatial amplitude window (wave packet)
packet_center_x = -400e3;
packet_center_y = -400e3;
packet_width_x = 4*400e3;
packet_width_y = 4*300e3;

% Synthetic wave scaling factor for spatial coordinates
% This factor modifies DX, which represents the real-world distance per pixel.
% A smaller DX means that each pixel covers a smaller physical distance, 
% effectively "zooming in" and making the wave pattern appear larger in the image.
% Conversely, a larger DX means that each pixel covers a larger real-world distance, 
% making the wave pattern appear smaller and more compressed.
%
% Example: 
%   - DXFactor = 1 means the default scaling (1:1 with pixel size).
%   - DXFactor = 1/4 means the wave features are stretched, appearing 4x larger.
%   - DXFactor = 4 means the wave features are shrunk, appearing 4x smaller.
DXFactor = 1;  

% Seconds between frames (important for drift or wave stepping)
time_resolution = 1800;

%----------- VIDEO OUTPUT SETTINGS -----------------------
createOutputVideo = true;    % Set to true to generate a video of processed frames
videoFrameRate = 5;          % Frames per second for the output video
videoApplyShrink = false;     % Use true to apply shrinkfactor to video frames
videoApplyWindow = true;     % Use true to apply windowing to video frames

%----------- PEAK-DETECTION/Brightness decay SETTINGS -----------------------------------
speed_min_threshold   = 7.5;        % [m s-1] absolute floor
speed_std_factor      = 2;         % N·σ above local mean
maxPeaksPerROI        = 3;         % safety cap (set [] for unlimited)
clevfactor_real       = 1.5;       % divides contour levels (real)
clevfactor_imag       = 2;         % divides contour levels (imag)
contourOption       = 'percentile';           % 'percentile' | '3sigma'
contourArray        = [50 60 70];             % if 'percentile'

radStep_pixels      = 2;       % width of annuli (original-res pixels)
radMax_pixels       = 120;     % stop radius for B-decay

showQuicklooks      = true;    % one figure/ROI with masks & B(r)
saveQuicklooks      = true;   % if you want files instead of onscreen
%-----------------------------------------------------------------------

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

%% 2.1) OPTIONAL VIDEO GENERATION
if createOutputVideo && numFrames > 0 % Check if video creation is enabled and frames exist
    fprintf('\n--- Generating Output Video ---\n');

    singleOutDir = 'C:\Users\admin\Documents\GitHub\Stratocu-Waves-Jr\test';

    % Define output video filename (customize as needed)
    videoFileName = sprintf('%s_Video_%s_to_%s.mp4', ...
                        upper(instrument), ...
                        datestr(fTimes(1), 'yyyymmddHHMM'), ...
                        datestr(fTimes(end), 'yyyymmddHHMM'));
    outputVideoFile = fullfile(singleOutDir, videoFileName); % Save in the same test dir

    % Determine if shrinking/windowing should be applied based on script settings AND video flags
    applyShrinkVideo = videoApplyShrink && (shrinkfactor ~= 1);
    applyWindowVideo = videoApplyWindow && doWindow;

    % Call the video creation function
    createVideoFromFrames(fNames, fTimes, dataDir, varName, ...
                       outputVideoFile, videoFrameRate, ...
                       applyShrinkVideo, shrinkfactor, invshrinkfactor, ... % Pass shrink args
                       applyWindowVideo, windowType, radius_factor, decay_rate, ... % Pass window args
                       instrument, methodName, ... % Pass preprocessing args
                       IR_threshold, IR_fillPercentile, ...
                       VIS_lowerPercentile, VIS_upperPercentile, ...
                       VIS_fillPercentile, clipMinHP, clipMaxHP, ...
                       lowPassFilterWidth_20, lowPassFilterWidth_50, lowPassFilterWidth_100, ...
                       Insolation_Correction, ...
                       syntheticWaveMode, ...
                       driftMode, ...
                       time_resolution, ...
                       drift_speed_m_s, driftAngleDeg, ...
                       cphase, wavelength, direction, zamplitude, ...
                       PBLdepth, dB_dzPBL, ...
                       packet_center_x, packet_center_y, ...
                       packet_width_x, packet_width_y, ...
                       DXFactor);

elseif createOutputVideo
    fprintf('\nVideo creation skipped: No frames were processed.\n');
end

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

if shrinkfactor ~= 1
    pixel_size_km = original_px_km * shrinkfactor;
    Scales = Scales / shrinkfactor;
    scalesForAdvection = scalesForAdvection/shrinkfactor; 
else
    pixel_size_km = original_px_km;
end

for f_idx = 1:numFrames
    % ---- Prepare basic info about this frame ----
    thisTime = fTimes(f_idx);
    frameDateStr = datestr(thisTime, 'yyyy_mm_dd_HHMMSS');
    fprintf('\nProcessing frame [%d/%d]: %s\n', f_idx, numFrames, frameDateStr);

    %singleOutDir = fullfile(rootSepacDir, 'Test');
    singleOutDir = 'C:\Users\admin\Documents\GitHub\Stratocu-Waves-Jr\test';
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
            DX = 1000 * original_px_km * DXFactor;
            Xm = (X - mean(X(:))) * DX;
            Ym = (Y - mean(Y(:))) * DX * -1;  % negative if Y runs downward
        end

    elseif syntheticWaveMode || driftMode
        % If doing synthetic wave or drift, start from the same base frame each iteration:
        data = base_frame;
    end

    % ---- Synthetic wave injection (if enabled) ----
    if syntheticWaveMode

        % --- meteo → math angle conversion (0° = Est, trigonometrical direction) ---
        theta = deg2rad(90 - direction);    % direction given in meteorological convention
        
        % --- wave vector ---
        k = (2*pi / wavelength) * cos(theta);   % kx
        l = (2*pi / wavelength) * sin(theta);   % ky  (kept for phase)
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


        % --- horizontal displacement associayted to w' ---
        dxy = (zamplitude / PBLdepth) * wavelength .* ...
               sin(phase - pi/2) ./ DX;          % scalar amplitude in px
        dx  = dxy .* cos(theta).* Ampwindow;                 % composante O-E  (columns)
        dy  = dxy .* sin(theta).* Ampwindow;                 % composante S-N  (lines)
        
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
        pxShift = driftDistance_km / (original_px_km);  % e.g. if pixel_size_km=9 => pxShift=2000/9
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

    % ---- Preprocessing (thresholds, highpass, etc.) ----
    data_pre = preprocessFrame(data, instrument, methodName, ...
                thisFullPath, thisTime, ...
                IR_threshold, IR_fillPercentile, ...
                VIS_lowerPercentile, VIS_upperPercentile, ...
                VIS_fillPercentile, clipMinHP, clipMaxHP, ...
                lowPassFilterWidth_20, lowPassFilterWidth_50, lowPassFilterWidth_100, ...
                Insolation_Correction);

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

end  

%% 4A) ROI-BASED SUMMARIES

% Build ROI squares.

effective_degrees_per_pixel = degrees_per_pixel * shrinkfactor;
square_size_px = round(square_size_deg / effective_degrees_per_pixel);

x_buffer_range = (window_buffer+1) : (colsF - window_buffer);
y_buffer_range = (window_buffer+1) : (rowsF - window_buffer);

adjusted_frame_width  = length(x_buffer_range);
adjusted_frame_height = length(y_buffer_range);

num_squares_x = ceil(adjusted_frame_width  / square_size_px);
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
                        Scales_orig, Angles, DisplayValuePower, data_bg_pre, singleOutDir, 'Global PowerWaveRose Overlay.png');

produceOverlayWaveRose('Coherence', roiCohCell, squares, num_squares_x, num_squares_y, ...
                        Scales_orig, Angles, DisplayValueCoherence, data_bg_pre, singleOutDir, 'Global CoherenceWaveRose Overlay.png');

produceOverlayWaveRose('Phase', roiSpeedCell, squares, num_squares_x, num_squares_y, ...
                        Scales_orig, Angles, DisplayValuePhase, data_bg_pre, singleOutDir, 'Global PhaseWaveRose Overlay.png');

produceOverlayWaveRose('Speed', roiSpeedCell, squares, num_squares_x, num_squares_y, ...
                        Scales_orig, Angles, DisplayValueSpeed, data_bg_pre, singleOutDir, 'Global SpeedWaveRose Overlay.png');


%% 4C) Final Correction of Large-Scale Speeds

% Apply the correction on the ROI-level wave–rose speeds
roiSpeedCell_corrected = limitSpeedByScale(roiSpeedCell, Scales_orig, nyquistScales, upperCutoff, beta, decayFactor, decaySharpness, matrixMode);

% Then, you can use produceOverlayWaveRose to visualize the final results:
produceOverlayWaveRose('Speed Corrected', roiSpeedCell_corrected, squares, num_squares_x, num_squares_y, ...
    Scales_orig, Angles, DisplayValueSpeed, data_bg_pre, singleOutDir, 'Global Corrected SpeedWaveRose Overlay.png');

%% 4D) ADVECTION CORRECTION ON ROI-BASED ROSES

if doAdvectionEstimation
    [roiSpeedCell_corrected, roiResidCell] = ...
        applyAdvectionCorrectionROI(roiSpeedCell, roiCohCell, ...
                                   Scales, Angles, ...
                                   scalesForAdvection, cohThreshold, ...
                                   pixel_size_km, time_resolution);

    % Apply the correction on the ROI-level wave–rose speeds (for instance, on the advection-corrected data)
roiSpeedCell_corrected = limitSpeedByScale(roiSpeedCell_corrected, Scales, nyquistScales, upperCutoff, beta, decayFactor, decaySharpness, matrixMode);


produceOverlayWaveRose('Speed UnAdvected', roiSpeedCell_corrected, squares, num_squares_x, num_squares_y, ...
                        Scales_orig, Angles, DisplayValueSpeed, data_bg_pre, singleOutDir, 'Global UnAdvected SpeedWaveRose Overlay.png');

end

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
produceAggregatedWaveRose('Power Global', globalPower, Scales_orig, Angles, ...
    singleOutDir, 'Global PowerWaveRose', saverose, DisplayValuePower);

% ----- (B) Global Coherence Wave–Rose -----
globalCross = squeeze( sum(sum(crossSpec_sum, 1, 'omitnan'), 2, 'omitnan') );
globalB1    = squeeze( sum(sum(B1_auto_sum,   1, 'omitnan'), 2, 'omitnan') );
globalB2    = squeeze( sum(sum(B2_auto_sum,   1, 'omitnan'), 2, 'omitnan') );
% Compute coherence (using standard formula: |S12|^2 / (S11*S22)):
globalCoherence = ( abs(globalCross).^2 ) ./ (globalB1 .* globalB2 );
produceAggregatedWaveRose('Coherence Global', globalCoherence, Scales_orig, Angles, ...
    singleOutDir, 'Global CoherenceWaveRose', saverose, DisplayValueCoherence);

% ----- (C) Global Speed (Phase) Wave–Rose -----
% For phase differences, we average the complex exponentials (for circular averaging)
globalPhaseExp = squeeze( sum(sum(phaseExp_sum, 1, 'omitnan'), 2, 'omitnan') );
% Average over the number of pairs and spatial domain:
globalAvgPhase = angle( globalPhaseExp / (numPairs * numPixels) );
produceAggregatedWaveRose('Speed Global', globalAvgPhase, Scales_orig, Angles, ...
    singleOutDir, 'Global SpeedWaveRose', saverose, DisplayValueSpeed);

% ----- (D) Global Phase Wave–Rose -----
% For phase differences, we average the complex exponentials (for circular averaging)
produceAggregatedWaveRose('Phase Global', globalAvgPhase, Scales_orig, Angles, ...
    singleOutDir, 'Global PhaseWaveRose', saverose, DisplayValuePhase);

% ----- (E) Global Corrected Speed (Phase) Wave–Rose -----
% Apply the correction on the global level wave–rose speeds
globalAvgPhase_corrected = limitSpeedByScale({globalAvgPhase}, Scales_orig, nyquistScales, upperCutoff,  beta, decayFactor, decaySharpness, matrixMode);
produceAggregatedWaveRose('Speed Corrected Global', cell2mat(globalAvgPhase_corrected), Scales_orig, Angles, ...
    singleOutDir, 'Global Corrected SpeedWaveRose', saverose, DisplayValueSpeed);

% ----- (F) Global Unadvected Speed (Phase) + Correction Wave–Rose -----
if doAdvectionEstimation
    [globalAvgPhase_unadvected, ResidCell] = ...
        applyAdvectionCorrectionROI(globalAvgPhase_corrected, {globalCoherence}, ...
                                   Scales, Angles, ...
                                   scalesForAdvection, cohThreshold, ...
                                   pixel_size_km, time_resolution);

    produceAggregatedWaveRose('Speed UnAdvected', cell2mat(globalAvgPhase_unadvected),Scales_orig, Angles, ...
        singleOutDir,  'Global UnAdvected SpeedWaveRose.png', saverose, DisplayValueSpeed);

end

fprintf('\nAll done. Single–frame and cross–temporal wavelet (coherence) computations complete.\n');

%% 6) PEAK-OVERLAY – draw real & imaginary contours of gravity-wave peaks
% ----------BEFORE STARTUP, BUILD ONE ANNOTATION CUBE FOR THE WHOLE FRAME ----------
Ny = size(data_bg_pre,1);   % full-res dims (after any resize/window)
Nx = size(data_bg_pre,2);

crestMask2D  = false(Ny, Nx); % 1 = crest
troughMask2D = false(Ny, Nx); % 1 = trough
realSum2D     = zeros (Ny, Nx, 'single'); % signed Σ(real(coeff_full))

spec_vis = spec_full;          % size = Ny × Nx × NSCALES × NANGLES

% Make an empty RGB image that we can draw on
figPeak = figure('visible','on');
imagesc(imresize(data_bg_pre, invshrinkfactor)); axis image off; colormap(gray); hold on;

% Loop over ROIs
for iROI = 1:numSquares
    % Convenience aliases
    speedMat = roiSpeedCell_corrected{iROI};          % [NSCALES × NANGLES]
    % For each row (corresponding to a fine-scale value), multiply by that scale.
    for iRow = 1:size(speedMat,1)
         speedMat(iRow, :) = (pixel_size_km*1000*(speedMat(iRow, :)/(2*pi)) .* Scales(iRow)* pi/sqrt(2))/(1800*shrinkfactor);
    end
    xR       = squares(iROI).x_range;
    yR       = squares(iROI).y_range;

    % Build a logical mask of peak candidates --------------------------
    mu   = mean(speedMat(:),'omitnan');
    sigma= std( speedMat(:),'omitnan');
    thr  = max(speed_min_threshold, mu + speed_std_factor*sigma);
    
    speedVals = abs(speedMat);

    visited   = false(size(speedVals));
    
    % we will collect blobs in a struct array
    blobList  = struct('members',{},'maxAmp',{});   

    while true
        % 1) find the strongest un-visited seed -----------------------------
        speedVals(visited) = -inf;                        % mask out visited
        [peakAmp, linearIdx] = max(speedVals(:));
        if peakAmp < thr, break, end                    % nothing above threshold
    
        [sSeed,aSeed] = ind2sub(size(speedVals), linearIdx);
    
        % 2) flood-fill ------------------------------------------------------
        thisBlob  = [];               % list of (s,a)
        queue     = [sSeed, aSeed];
    
        while ~isempty(queue)
            s = queue(1,1);  a = queue(1,2);
            queue(1,:) = [];                         % pop
            if visited(s,a), continue, end
            visited(s,a) = true;                     % mark now
    
            if abs(speedMat(s,a)) < thr, continue, end
    
            thisBlob(end+1,:) = [s,a];               %#ok<AGROW>
    
            % enqueue 8-neighbours with angle wrap
            for ds = -1:1
                for da = -1:1
                    if ds==0 && da==0, continue, end
                    ss = s + ds;
                    aa = mod(a - 1 + da, NANGLES) + 1;
                    if ss>=1 && ss<=NSCALES && ~visited(ss,aa)
                        queue(end+1,:) = [ss, aa];    %#ok<AGROW>
                    end
                end
            end
        end
    
        % 3) aggregate complex coefficients over this blob ------------------
        wav_c_blob = zeros(round(Ny/shrinkfactor), round(Nx/shrinkfactor), 'like', spec_vis(:,:,1,1));
        for jj = 1:size(thisBlob,1)
            sIdx = thisBlob(jj,1);   aIdx = thisBlob(jj,2);
            wav_c_blob = wav_c_blob + spec_vis(:,:,sIdx,aIdx);
        end
        wav_real = real(wav_c_blob);
        wav_imag = imag(wav_c_blob);

        % Optionally restrict the display to the current ROI only ----
        mask = false(size(wav_real));
        mask(yR, xR) = true;
        wav_real(~mask) = NaN;
        wav_imag(~mask) = NaN;

        % Choose contour levels based on either absolute values or percentiles.
        switch lower(contourOption)
            case 'percentile'
                contourLevels = prctile(abs(wav_real(:)), contourArray);
                % Plot contours:
                % For positive values:
                contour(wav_real, contourLevels, 'LineColor', 'r', 'LineWidth', 0.5);
                % For negative values (mirror the levels):
                contour(wav_real, -contourLevels, 'LineColor', 'b', 'LineWidth', 0.5);
            case '3sigma'
                % Robust contour levels based on local std --------------------
                sig_r = std(wav_real(:),'omitnan');
                sig_i = std(wav_imag(:),'omitnan');
        
                posLevels =  (sig_r/clevfactor_real) : (sig_r/clevfactor_real) :  3*sig_r;
                negLevels = -(sig_r/clevfactor_real) :-(sig_r/clevfactor_real) : -3*sig_r;
                imagLevels=  (sig_i/clevfactor_imag) : (sig_i/clevfactor_imag) :  3*sig_i;
        
                %Plot: red = +real, blue = –real, green dashed = |imag| ------
                contour(wav_real, posLevels, 'LineColor','red',  'LineWidth',0.7);
                contour(wav_real, negLevels, 'LineColor','blue', 'LineWidth',0.7);
                %contour(wav_imag, imagLevels,'LineColor',[0 0.6 0], ...
                %                               'LineStyle','--','LineWidth',0.7);

            otherwise
                error('Unknown contour option. Choose either "3sigma" or "percentile".');
        end

        %=====  up-sample to full resolution  ======================
        if shrinkfactor ~= 1
            coeff_full  = imresize(wav_c_blob,   shrinkfactor, 'bilinear');  % complex
            mask_full   = imresize(mask,    shrinkfactor, 'nearest')>0;
            coeff_full = coeff_full(1:end-1,:);
            mask_full = mask_full(1:end-1,:);
        else
            coeff_full  = wav_c_blob;   % already full-res
            mask_full   = mask;
        end
        
        % Recaculate the level like above
        switch lower(contourOption)
          case 'percentile'
            lev = prctile(abs(real(coeff_full(mask_full))), contourArray(1));
          case '3sigma'
            lev = std(real(coeff_full(mask_full)), 'omitnan');
        end
        
        % create the two 2D masks
        crest  = ( real(coeff_full) >=  lev ) & mask_full;
        trough = ( real(coeff_full) <= -lev ) & mask_full;
        
        % Merges with the global masks
        crestMask2D  = crestMask2D  | crest;
        troughMask2D = troughMask2D | trough;
        realSum2D     = realSum2D  +  single(real(coeff_full) .* mask_full);                   
       
        % optional: keep a list of blobs for later ranking
        blobList(end+1).members = thisBlob;            %#ok<AGROW>
        blobList(end).maxAmp   = peakAmp;
    end

    % Draw ROI rectangle so the user knows which square is which -------
    rectangle('Position',[xR(1), yR(1), ...
               numel(xR), numel(yR)], ...
               'EdgeColor','k','LineWidth',0.8);
end

% Final cosmetics & export ---------------------------------------------
title(sprintf('Peak wavelet contours – %s', datestr(thisTime)), ...
      'Color','w','FontWeight','bold');

peakFileName = fullfile(singleOutDir, sprintf('Peaks_%s.png', ...
                     datestr(thisTime,'yyyymmdd_HHMMSS')));
exportgraphics(figPeak, peakFileName, 'Resolution', 300);
close(figPeak);

fprintf('Saved peak-overlay to: %s\n', peakFileName);

% --- NetCDF export --------------------------------
fname = sprintf('%s_spectralpeaks_%s_to_%s.nc', ...
                upper(instrument), ...
                datestr(startDate,'yyyymmddHHMM'), ...
                datestr(endDate,  'yyyymmddHHMM'));
outNC = fullfile(singleOutDir, fname);

if exist(outNC,'file'), delete(outNC); end

nccreate(outNC,'REAL_SUM', 'Datatype','single', ...
                     'Dimensions',{'y',Ny,'x',Nx});
ncwrite (outNC,'REAL_SUM', single(realSum2D));

ncwriteatt(outNC,'/','instrument' ,instrument);

ncwriteatt(outNC,'/','time_coverage_start', char(startDate));
ncwriteatt(outNC,'/','time_coverage_end',   char(endDate));

ncwriteatt(outNC,'/','n_scales', int32(NSCALES));
ncwriteatt(outNC,'/','n_angles', int32(NANGLES));

ncwriteatt(outNC,'/','description', ...
   'Accumulated real part of selected CWT peaks');
fprintf('   ↳ wrote REAL_SUM → %s\n', outNC);



%% 6.1) BRIGHTNESS–DECAY ANALYSIS  ________________________________________
%
%  • detects each spectral blob (flood-fill in scale/angle space)
%  • builds a ridge (skeleton) along the positive real-part lobe
%  • chooses an anchor pixel (brightest on ridge) – for thumbnails only
%  • computes B(r) = ⟨ΔT⟩(r) as a function of distance to the whole ridge
%  • stacks all B(r) profiles → composite curve
%  • optional: postage-stamp stack around every anchor
%  ________________________________________________________________________


spec_vis   = spec_full;                     % Ny×Nx×NSCALES×NANGLES
data_full  = data_bg_pre;                   % full-res, pre-processed
[Ny_full, Nx_full] = size(data_full);

% radial bins for B(r)
radBins   = (0:radStep_pixels:radMax_pixels).';     % column vector
nBins     = numel(radBins) - 1;
Bstack    = [];                                     % Npeaks × nBins

% options for § 6.2
doPostageStamp   = true;
doAzimAvg        = true;
stampHalfSize_px = radMax_pixels;
saveStampNetCDF  = true;

StackStamp  = [];          % (2H+1)×(2H+1)×N
StampCtrXY  = [];          % N×2 (cx,cy)
StampRadProf= [];          % N×nBins

NSCALES  = numel(Scales);
NANGLES  = numel(Angles);

% -------------------------------------------------------------------------
%                           LOOP OVER ROIs
% -------------------------------------------------------------------------
for iROI = 1:numSquares

    % ---------- A) speed matrix and threshold ---------------------------
    speedMat = roiSpeedCell_corrected{iROI};
    for iRow = 1:size(speedMat,1)
        speedMat(iRow,:) = ...
            (pixel_size_km*1000*(speedMat(iRow,:)/(2*pi)) .* ...
             Scales(iRow) * pi/sqrt(2)) / (1800*shrinkfactor);
    end

    xR_sh = squares(iROI).x_range;
    yR_sh = squares(iROI).y_range;

    mu  = mean(speedMat(:),'omitnan');
    sg  = std (speedMat(:),'omitnan');
    thr = max(speed_min_threshold, mu + speed_std_factor*sg);

    speedAbs = abs(speedMat);
    visited  = false(size(speedAbs));

    % ROI indices in full resolution
    xR_full = ((xR_sh-1)*shrinkfactor)+1 : (xR_sh(end)*shrinkfactor);
    yR_full = ((yR_sh-1)*shrinkfactor)+1 : (yR_sh(end)*shrinkfactor);

    % ---------------------------------------------------------------------
    %                     LOOP OVER BLOBS  (one blob = one crest)
    % ---------------------------------------------------------------------
    while true
        speedAbs(visited) = -inf;
        [peakAmp, linIdx] = max(speedAbs(:));
        if peakAmp < thr, break, end

        % ---------- quick-look canvas ---------------------------------------
        if showQuicklooks
            figQL = figure('Name',sprintf('ROI %d',iROI), ...
                           'Units','normalized','Position',[.05 .05 .35 .6]);
            subplot(2,2,1);
            imagesc(imresize(data_full,invshrinkfactor));
            axis image off; colormap gray; hold on
            rectangle('Position',[xR_sh(1) yR_sh(1) ...
                                  numel(xR_sh) numel(yR_sh)], ...
                      'EdgeColor','c','LineWidth',1.5);
            title(sprintf('ROI %d context',iROI));
        end

        [sSeed,aSeed] = ind2sub(size(speedAbs), linIdx);

        % --- flood-fill in (scale,angle) space --------------------------
        thisBlob = [];                        % list of (s,a)
        queue    = [sSeed, aSeed];
        while ~isempty(queue)
            s = queue(1,1); a = queue(1,2); queue(1,:) = [];
            if visited(s,a), continue, end
            visited(s,a) = true;
            if abs(speedMat(s,a)) < thr, continue, end

            thisBlob(end+1,:) = [s,a];                  %#ok<AGROW>

            for ds = -1:1
                for da = -1:1
                    if ds==0 && da==0, continue, end
                    ss = s + ds;
                    aa = mod(a - 1 + da, NANGLES) + 1;
                    if ss>=1 && ss<=NSCALES && ~visited(ss,aa)
                        queue(end+1,:) = [ss, aa];      %#ok<AGROW>
                    end
                end
            end
        end

        % --- aggregate complex coefficients over the blob --------------
        wav_c_blob = zeros(size(spec_vis(:,:,1,1)),'like',spec_vis);
        for jj = 1:size(thisBlob,1)
            wav_c_blob = wav_c_blob + ...
                         spec_vis(:,:, thisBlob(jj,1), thisBlob(jj,2));
        end
        wav_r = real(wav_c_blob);

        % restrict to current ROI
        maskROI = false(size(wav_r));
        maskROI(yR_sh,xR_sh) = true;
        wav_r(~maskROI) = NaN;

        % positive-lobe threshold in spatial domain
        switch lower(contourOption)
            case 'percentile'
                posLevel = prctile(abs(wav_r(:)), contourArray(1));
            case '3sigma'
                posLevel = std(wav_r(:),'omitnan');
            otherwise
                error('Unknown contour option');
        end
        maskRed = wav_r >= posLevel;
        maskRed(~maskROI) = false;

        % keep largest spatial blob
        CCred = bwconncomp(maskRed,4);
        if CCred.NumObjects==0, continue, end
        [~,bigI] = max(cellfun(@numel,CCred.PixelIdxList));
        maskRed  = false(size(maskRed));
        maskRed(CCred.PixelIdxList{bigI}) = true;

        % --- ridge (skeleton) extraction --------------------------------
        ridgeMask = bwskel(maskRed);
        if ~any(ridgeMask(:)), continue, end
        [ySkel,xSkel] = find(ridgeMask);

        % map ridge to full resolution
        ridgeMask_full = false(Ny_full,Nx_full);
        row_full = min(max(round(ySkel*shrinkfactor),1),Ny_full);
        col_full = min(max(round(xSkel*shrinkfactor),1),Nx_full);
        ridgeMask_full(sub2ind([Ny_full Nx_full],row_full,col_full)) = true;

        % anchor = ridge pixel with max ΔT
        subROI_full = data_full(yR_full,xR_full);
        baselineROI = mean(subROI_full(:),'omitnan');
        deltaT_full = data_full - baselineROI;
        deltaT_on_ridge = deltaT_full(ridgeMask_full);
        [~,bestIdx] = max(deltaT_on_ridge);
        cx = col_full(bestIdx);   cy = row_full(bestIdx);


        % --- ridge-centred brightness-decay profile ---------------------
        Rmap = bwdist(ridgeMask_full);          % distance [px] to *entire* ridge
        
        prof = nan(1,nBins);
        for jj = 1:nBins
            annulus = Rmap>=radBins(jj) & Rmap<radBins(jj+1);
            if any(annulus(:))
                prof(jj) = mean(deltaT_full(annulus),'omitnan');
            end
        end
        if all(isnan(prof)), continue, end      % <-- keep guard
        Bstack = [Bstack; prof];                %#ok<AGROW>


        % --- optional postage-stamp ------------------------------------
        if doPostageStamp
            H = stampHalfSize_px;  S = 2*H+1;
            rowL = max(cy-H,1); rowR = min(cy+H,Ny_full);
            colL = max(cx-H,1); colR = min(cx+H,Nx_full);

            stampSq  = NaN(S,S,'like',deltaT_full);
            rowDst   = 1 + (rowL-(cy-H));
            colDst   = 1 + (colL-(cx-H));
            stampSq(rowDst:rowDst+(rowR-rowL), ...
                    colDst:colDst+(colR-colL)) = ...
                    deltaT_full(rowL:rowR,colL:colR);

            StackStamp(:,:,end+1) = stampSq;            %#ok<AGROW>
            StampCtrXY(end+1,:)   = [cx cy];            %#ok<AGROW>

            if doAzimAvg
                prof_stamp = azimProfile(stampSq,H+1,H+1,radBins);
                if all(isnan(prof_stamp)),  continue,  end    % <-- guard
                StampRadProf = [StampRadProf; prof_stamp]; %#ok<AGROW>
            end
        end

        % --- quick-look panels -----------------------------------------
        if showQuicklooks
            subplot(2,2,2); hold on
            axis ij; 
            imagesc(imresize(data_full,invshrinkfactor)); axis image off
            contour(wav_r,[posLevel posLevel],'r','LineWidth',.7);
            [pY,pX] = find(ridgeMask);
            plot(pX,pY,'r.', 'MarkerSize',4);
            xlim([xR_sh(1) xR_sh(end)]);
            ylim([yR_sh(1) yR_sh(end)]);
            title('CWT real + ridge');

            subplot(2,2,3);
            imagesc(data_full); axis image off; colormap gray
            xlim([xR_full(1) xR_full(end)]);
            ylim([yR_full(1) yR_full(end)]);
            title('Full-res & anchor');

            subplot(2,2,4);
            rPlot = (radBins(1:end-1)+radBins(2:end))/2;
            plot(rPlot,prof,'-o'); grid on
            xlabel('r [px]'); ylabel('\DeltaT');
            title('Brightness-decay');
        end

        % save or show quick-look
        if showQuicklooks
            if saveQuicklooks
                qlName = fullfile(singleOutDir, ...
                         sprintf('ROI%02d_quicklook_%s.png', iROI, ...
                         datestr(thisTime,'yyyymmdd_HHMMSS')));
                exportgraphics(figQL, qlName, 'Resolution',250);
                close(figQL)
            else
                drawnow
            end
        end

    end % while blobs

end % iROI

% ---------- 6.1.1) COMPOSITE BRIGHTNESS-DECAY ---------------------------
if ~isempty(Bstack)
    Bcomp = nanmean(Bstack,1);
    figure('Name','Composite brightness-decay', ...
           'Units','normalized','Position',[.45 .3 .25 .4]);
    rPlot = (radBins(1:end-1)+radBins(2:end))/2;
    plot(rPlot,Bcomp,'LineWidth',2); grid on
    xlabel('Radius r [original pixels]');
    ylabel('\langle\DeltaT\rangle');
    title(sprintf('Composite B-decay over %d lobes', size(Bstack,1)));
end


%% 6.2  POSTAGE-STAMP STACKING ________________________

if doPostageStamp && ~isempty(StackStamp)
    H = stampHalfSize_px;        % shorthand
    N = size(StackStamp,3);

    StampMean = nanmean(StackStamp,3);
    StampStd  =  nanstd(StackStamp,0,3);

    % -- quick-look --
    figure('Name','Postage-stamp composite')
    subplot(1,2,1); imagesc(StampMean); axis image off
    title(sprintf('\\DeltaT mean of %d stamps',N));
    colormap(parula)
    subplot(1,2,2); imagesc(StampStd); axis image off
    title(sprintf('\\DeltaT std-dev'));
    drawnow

    % -- optional composite radial profile (B2) --
    if doAzimAvg && ~isempty(StampRadProf)
        figure('Name','Composite B(r) from stamps')
        rPlot = (radBins(1:end-1)+radBins(2:end))/2;
        plot(rPlot, nanmean(StampRadProf,1),'LineWidth',2); grid on
        xlabel('r [px]'); ylabel('\langle\DeltaT\rangle');
        title(sprintf('Mean of %d stamps',N));
    end

    % -- NetCDF export ------------------------------------------
    if saveStampNetCDF
        fname = sprintf('%s_brightsnaps_%s_to_%s.nc', ...
                        upper(instrument), ...
                        datestr(startDate,'yyyymmddHHMM'), ...
                        datestr(endDate,  'yyyymmddHHMM'));
        outNC = fullfile(singleOutDir, fname);
        if exist(outNC,'file'); delete(outNC); end

        nccreate(outNC,'STAMP_MEAN','Datatype','single', ...
                           'Dimensions',{'y',2*H+1,'x',2*H+1});
        ncwrite (outNC,'STAMP_MEAN',single(StampMean));

        nccreate(outNC,'STAMP_STD', 'Datatype','single', ...
                           'Dimensions',{'y',2*H+1,'x',2*H+1});
        ncwrite (outNC,'STAMP_STD', single(StampStd));

        % optional – save the stack itself (big!)
        % nccreate(outNC,'STAMP_STACK','Datatype','single', ...
        %     'Dimensions',{'y',2*H+1,'x',2*H+1,'k',N});
        % ncwrite (outNC,'STAMP_STACK',single(StackStamp));

        ncwriteatt(outNC,'/','description', ...
            'Postage-stamp stack around Brightness Points - mean & std');

        ncwriteatt(outNC,'/','time_coverage_start', char(startDate));
        ncwriteatt(outNC,'/','time_coverage_end',   char(endDate));
        
        ncwriteatt(outNC,'/','n_stamps',int32(N));
        fprintf('   ↳ wrote POSTAGE-STAMP composites → %s\n',outNC);
    end
end

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
    if evalin('base','shrinkfactor')==2
        sizeFactor = 2;    % Factor to enlarge each inset.
        xOffset   = 0.035; % Horizontal offset between insets.
        yOffset   = 0.08;  % Vertical offset between insets.
    elseif evalin('base','shrinkfactor')==1
        sizeFactor = 1;
        xOffset   = -0.027; % Horizontal offset between insets.
        yOffset   = 0;  % Vertical offset between insets.
    end
    
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
            if evalin('base','shrinkfactor')==2
                xNorm = posMain(1) + ((xMin - xLim(1)) / (xLim(2) - xLim(1))) * posMain(3);
                wNorm = (width / (xLim(2) - xLim(1))) * posMain(3) * sizeFactor;
                yNorm = posMain(2) + ((yLim(2) - yMax) / (yLim(2) - yLim(1))) * posMain(4);
                hNorm = (height / (yLim(2) - yLim(1))) * posMain(4) * sizeFactor;
            elseif evalin('base','shrinkfactor')==1
                % Convert ROI data coordinates to figure–normalized coordinates.
                xNorm = posMain(1) + ((xMin - xLim(1)) / (xLim(2) - xLim(1))) * posMain(3)*sizeFactor;
                wNorm = (width / (xLim(2) - xLim(1))) * posMain(3) * sizeFactor;
                yNorm = posMain(2) + ((yLim(2) - yMax) / (yLim(2) - yLim(1))) * posMain(4)*sizeFactor;
                hNorm = (height / (yLim(2) - yLim(1))) * posMain(4) * sizeFactor;
            end
            
            % Adjust for spacing between insets.
            xNorm = xNorm + (ix - 1) * xOffset;
            yNorm = yNorm - (iy - 1) * yOffset;
            
            % Create an inset axes at the computed position.
            if evalin('base','shrinkfactor')==2
                axInset = axes('Position', [xNorm+0.07, yNorm-0.08, wNorm, hNorm]);
            elseif evalin('base','shrinkfactor')==1
                axInset = axes('Position', [xNorm+0.07, yNorm, wNorm, hNorm]);
            end

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
    
    % --- Compute edges for angles and log scales ---
    % Angular edges (midpoints between given angles)
    d_angle = diff(Angles(1:2)); % Assumes uniform spacing
    angle_edges = [Angles - d_angle/2, Angles(end) + d_angle/2];
    %angle_edges = Angles - d_angle/2;

    % Log-scale edges (midpoints between log10(Scales))
    log_Scales = log10(Scales);
    d_log = diff(log_Scales(1:2)); % Assumes uniform log spacing
    log_edges = [log_Scales(1) - d_log/2, log_Scales + d_log/2];
    
    % Create meshgrid of edges
    [Theta_edges, R_edges] = meshgrid(angle_edges, log_edges);

   
    if startsWith(labelStr, 'Speed', 'IgnoreCase', true)
        pixel_size_km = evalin('base','pixel_size_km');
        shrinkfactor = evalin('base','shrinkfactor');
        % For each row (corresponding to a fine-scale value), multiply by that scale.
        for iRow = 1:size(innerpower,1)
             innerpower(iRow, :) = (pixel_size_km*1000*(innerpower(iRow, :)/(2*pi)) .* Scales(iRow)* pi/sqrt(2))/(1800*shrinkfactor);
        end
    end
    
    % Convert polar coordinates to Cartesian (using a log-scale for R).
    [X_pos, Y_pos] = pol2cart(Theta_edges, R_edges);
    padded_power = padarray(innerpower, [1 1], NaN, 'post'); % Add NaNs to edges
   
    % Plot the positive half of the rose.
    pcolor(axHandle,X_pos, Y_pos, padded_power);
    shading(axHandle, 'flat');
    
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
                val_linear = 10^(level_log);
        %text(axHandle, level_log, 0, sprintf('%.2f', val_linear), 'Color','k','FontSize',5,'HorizontalAlignment','left');
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
        maxVal = [];  % Auto-scale if not provided
    end

    if size(maxVal,1)==2
        minVal = maxVal(1);
        maxVal = maxVal(2);
    end

    innerpower = waveRoseMat;

    % --- Compute edges for angles and log scales ---
    % Angular edges (midpoints between given angles)
    d_angle = diff(Angles(1:2)); % Assumes uniform spacing
    angle_edges = [Angles - d_angle/2, Angles(end) + d_angle/2];
    %angle_edges = Angles - d_angle/2;

    % Log-scale edges (midpoints between log10(Scales))
    log_Scales = log10(Scales);
    d_log = diff(log_Scales(1:2)); % Assumes uniform log spacing
    log_edges = [log_Scales(1) - d_log/2, log_Scales + d_log/2];
    
    % Create meshgrid of edges
    [Theta_edges, R_edges] = meshgrid(angle_edges, log_edges);

    % --- Modification for Speed ---
    if startsWith(labelStr, 'Speed', 'IgnoreCase', true)
        pixel_size_km = evalin('base','pixel_size_km');
        shrinkfactor = evalin('base','shrinkfactor');
        % For each row (corresponding to a fine-scale value), multiply by that scale.
        for iRow = 1:size(innerpower,1)
             innerpower(iRow, :) = (pixel_size_km*1000*(innerpower(iRow, :)/(2*pi)) .* Scales(iRow)* pi/sqrt(2))/(1800*shrinkfactor);
        end
    end

    % Convert to Cartesian coordinates
    [X_pos, Y_pos] = pol2cart(Theta_edges, R_edges);
    
    padded_power = padarray(innerpower, [1 1], NaN, 'post'); % Add NaNs to edges
   
    % Create figure and plot
    figRose = figure('visible','off');
    set(figRose, 'Position', [100 100 600 600]);
    ax1 = axes('Position',[0.1 0.1 0.75 0.75]);
    hold(ax1, 'on');
    pcolor(ax1, X_pos, Y_pos, padded_power);
    shading(ax1, 'flat');
    
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
        clim(ax1, [0, max(innerpower(:))]);
    end

    % Set hard axis limits to contain upper half
    xlim(ax1, [min(X_pos,[],'all')*1.1  max(X_pos,[],'all')*1.1 ]);
    ylim(ax1, [min(Y_pos,[],'all')*1.1 max(Y_pos,[],'all')*1.05]);  % Restrict to upper half

    axis(ax1, 'equal', 'tight', 'off');

    % Drawing outlines
    theta_outline = linspace(0-(d_angle/2), pi+(d_angle/2), 100);  
    [x_outline, y_outline] = pol2cart(theta_outline, max(X_pos,[],'all'));
    plot(ax1, x_outline, y_outline, 'k-', 'LineWidth', 1.5);

    inner_radius = log_edges(1);
    theta_outline_inner = linspace(0 - d_angle/2, pi + d_angle/2, 100);
    [x_inner, y_inner] = pol2cart(theta_outline_inner, inner_radius);
    plot(ax1, x_inner, y_inner, 'k-', 'LineWidth', 1.5);

    [x_debut, y_debut] = pol2cart(pi + d_angle/2, inner_radius);
    [x_fin, y_fin] = pol2cart(pi + d_angle/2, log_edges(end));

    plot(ax1, [x_debut, x_fin], [y_debut, y_fin], 'k-', 'LineWidth', 1.5);

    [x_debut, y_debut] = pol2cart(0 - d_angle/2, inner_radius);
    [x_fin, y_fin] = pol2cart(0 - d_angle/2, log_edges(end));

    plot(ax1, [x_debut, x_fin], [y_debut, y_fin], 'k-', 'LineWidth', 1.5);
    
    % Add radial grid lines (using logarithmic scale for the scales).
    max_scale_log = log10(max(Scales));
    for i = 1:length(Scales)
        level_log = log10(Scales(i));
        theta_ring = linspace(0-(d_angle/2), pi+(d_angle/2), 180);
        [x_ring, y_ring] = pol2cart(theta_ring, level_log);
        plot(ax1, x_ring, y_ring, 'k--', 'LineWidth',0.5);
                val_linear = 10^(level_log);
        text(ax1, level_log, 0, sprintf('%.2f', val_linear), 'Color','k','FontSize',5,'HorizontalAlignment','left');
    end
    
    % Add angular lines.
    angle_ticks = linspace(0, pi, 7);
    max_r = max(X_pos,[],'all');
    for i = 1:length(angle_ticks)
        [x_label, y_label] = pol2cart(angle_ticks(i), max_r);
        line(ax1, [0 x_label], [0 y_label], 'Color',[0.5 0.5 0.5],'LineStyle','--');
    end
    
    title(ax1, sprintf('%s Wave–Rose', labelStr), 'FontSize',12, 'FontWeight','bold');
    
    c = colorbar(ax1, 'Location','eastoutside');
    c.Label.String = 'Wavelet Power';
    c.Label.FontWeight = 'bold';
    
    ax1_pos = ax1.Position;
    ax1_pos(3) = ax1_pos(3) * 0.85;
    ax1.Position = ax1_pos;
    
    if saverose
        roseFileName = fullfile(outDir, sprintf('%s.png', outName));
        fprintf('Saving wave–rose plot to: %s\n', roseFileName);
        exportgraphics(figRose, roseFileName, 'Resolution',300);
    end
    close(figRose);
end

%==========================================================================
function [u, alpha] = estimateAdvectionFromPhase(pdifrose, cohrose, Scales, Angles, ...
                                                 scalesForAdvection, cohThreshold, pixel_size_km, dt)
% estimateAdvectionFromPhase estimates the advection speed (u) and direction (alpha)
% from the averaged phase differences and coherence.
%
% Inputs:
%   pdifrose         : [NSCALES x NANGLES] matrix of mean phase differences (in radians)
%   cohrose          : [NSCALES x NANGLES] matrix of coherence (or coherence^2)
%   Scales           : vector of scales in pixels [NSCALES x 1]
%   Angles           : vector of angles in radians [1 x NANGLES]
%   scalesForAdvection: vector of scales to use (e.g. [4, 8, 16])
%   cohThreshold     : Coherence threshold; ignore bins with coherence < cohThreshold
%   pixel_size_km    : km per pixel (after scaling)
%   dt               : Time interval between frames (in seconds)
%
% Outputs:
%   u     : Estimated advection speed in m/s
%   alpha : Estimated advection direction (in radians)

% Select only the desired scales
mask_scales = ismember(Scales, scalesForAdvection);
pdif_sel = pdifrose(mask_scales, :);
coh_sel  = sqrt(cohrose(mask_scales, :));  % If cohrose is gamma^2, take the square root
S_sel    = Scales(mask_scales);

% Mask out bins where coherence is below the threshold
coh_mask = (coh_sel >= cohThreshold);
pdif_sel(~coh_mask) = NaN;

% Convert the scales (in pixels) into wavelengths in meters.
% (Assuming one cycle = π radians)
lambda_sel = S_sel * (pixel_size_km*1000) * (pi / sqrt(2));

% Use fminsearch to find the (u, alpha) that minimizes the phase error
x0 = [5, pi/2];  % initial guess: u=5 m/s, alpha=pi/2
opts = optimset('Display','none');
x_opt = fminsearch(@(x) costFun(x, pdif_sel, coh_sel, Angles, lambda_sel, dt), x0, opts);

u     = x_opt(1);
alpha = x_opt(2);
end

function err = costFun(x, pdif, coh, Angles, lambda, dt)
% costFun computes the weighted mean square error between observed phase and
% predicted phase for a given u and alpha.
%
% Inputs:
%   x     : [u, alpha]
%   pdif  : Observed phase differences [nSelScales x NANGLES]
%   coh   : Coherence weights for the selected bins [nSelScales x NANGLES]
%   Angles: Vector of angles (radians)
%   lambda: Vector of wavelengths corresponding to selected scales (in meters)
%   dt    : Time resolution (in seconds)
%
% Output:
%   err   : Mean squared error (weighted by coherence)

u     = x(1);
alpha = x(2);
[nSc, nAng] = size(pdif);
phi_pred = zeros(nSc, nAng);

% Compute the predicted phase for each scale and angle:
% phi_pred(s, theta) = ((u*dt)/lambda_s) * cos(theta - alpha) * pi
for iS = 1:nSc
    for iA = 1:nAng
        phi_pred(iS,iA) = ((u*dt)/lambda(iS)) * cos(Angles(iA)-alpha) * pi;
    end
end

% Compute the phase difference robustly using complex representation
diffPhase = angle(exp(1i*(pdif - phi_pred)));
% Weight by the coherence (or amplitude)
w = coh;
err = nanmean( (w(:).*diffPhase(:)).^2 );
end

function [roiSpeedCell_corrected, roiResidCell] = applyAdvectionCorrectionROI( ...
    roiSpeedCell, roiCohCell, ...
    Scales, Angles, ...
    scalesForAdvection, cohThreshold, ...
    pixel_size_km, dt)
% APPLYADVECTIONCORRECTIONROI
%   Loops over each ROI, retrieves the average phase (roiSpeedCell{iROI}),
%   estimates a global advection (u, alpha), then subtracts that advection's
%   predicted phase from the observed phase. Returns the corrected wave–rose
%   (roiSpeedCell_corrected) and an optional residual (roiResidCell).
%
%   Inputs:
%     roiSpeedCell  : cell array of [NSCALES x NANGLES] average phase
%     roiCohCell    : cell array of [NSCALES x NANGLES] coherence^2
%     Scales, Angles: wavelet scales (pixels) and angles (radians)
%     scalesForAdvection : subset of scales to use for the advection fit
%     cohThreshold       : ignore scale–angle bins with coherence < threshold
%     pixel_size_km      : real distance per pixel (km)
%     dt                 : time resolution (seconds)
%
%   Outputs:
%     roiSpeedCell_corrected : same size as roiSpeedCell, but with
%                              the large-scale advection subtracted out
%     roiResidCell           : an alternate "residual" measure (optional)

    numROI = numel(roiSpeedCell);
    roiSpeedCell_corrected = cell(size(roiSpeedCell));
    roiResidCell           = cell(size(roiSpeedCell));

    for iR = 1:numROI
        % 1) Extract the naive average phase (pdifrose) and coherence^2 for this ROI
        pdifrose = roiSpeedCell{iR};   % [NSCALES x NANGLES]
        cohrose  = roiCohCell{iR};     % [NSCALES x NANGLES], i.e. gamma^2

        % 2) Estimate the advection
        [u_hat, alpha_hat] = estimateAdvectionFromPhase( ...
            pdifrose, cohrose, ...
            Scales, Angles, ...
            scalesForAdvection, cohThreshold, ...
            pixel_size_km, dt);

        % 3) Build the theoretical advection phase for each scale, angle
        NSCALES  = numel(Scales);
        NANGLES  = numel(Angles);
        phi_adv  = zeros(NSCALES, NANGLES);
        for iS = 1:NSCALES
            lambda_m = Scales(iS) * (pixel_size_km*1000) * pi/sqrt(2);  % wavelength in meters
            for iA = 1:NANGLES
                phi_adv(iS,iA) = ((u_hat*dt)/lambda_m) * cos(Angles(iA) - alpha_hat) * pi;
            end
        end

        % 4) Subtract that from the observed phase (pdifrose)
        %    We'll do a robust difference: residualPhase = angle( e^{i( pdif - phi_adv )} ).
        diffPhase = pdifrose - phi_adv;
        residualPhase = angle( exp(1i * diffPhase) );

        % 5) Optionally store the "corrected" wave–rose
        %    You can interpret this as "phase minus advection". 
        roiSpeedCell_corrected{iR} = residualPhase;

        % 6) If desired, also store the difference in raw form
        %    (e.g., a simple difference in degrees)
        rawDiff = diffPhase;  % or convert to degrees, etc.
        roiResidCell{iR} = rawDiff; 
    end
end

%==========================================================================
function roiSpeedCell_final = limitSpeedByScale(roiSpeedCell_in, Scales, nyquistScales, upperCutoff, alpha, decayFactor, decaySharpness, matrixMode)
% LIMITSPEEDBYSCALE Corrects phase wrapping (Nyquist) and smooths large scales.
%
%   roiSpeedCell_final = limitSpeedByScale(roiSpeedCell_in, Scales, nyquistScales, upperCutoff, alpha)
%
% Inputs:
%   roiSpeedCell_in : Cell array of [NSCALES x NANGLES] phase matrices
%   Scales          : Vector of scale values
%   nyquistScales   : Scales to apply Nyquist correction (e.g., [1:4])
%   upperCutoff     : Upper scale cutoff for trend smoothing
%   alpha           : Smoothing factor for large scales [0,1]
%
% Output:
%   roiSpeedCell_final : Corrected phase data
% Version 2: Local Nyquist correction + constrained trend

nROI = numel(roiSpeedCell_in);
roiSpeedCell_final = cell(size(roiSpeedCell_in));
phaseJumpThreshold = pi * 0.6; % More conservative threshold
zeroPhaseThreshold = 1; % Threshold to detect near-zero phase points (radians)

for iR = 1:nROI
    currentWave = roiSpeedCell_in{iR};
    [nScales, nAngles] = size(currentWave);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 1: Local Nyquist Correction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nyquistMask = ismember(Scales, nyquistScales);    
    correctedWave = currentWave;

    % Calculate trusted profile from non-Nyquist scales
    trustedScalesMask = ~ismember(Scales, nyquistScales);
    trustedProfile = mean(currentWave(trustedScalesMask, :), 1);

    for s = find(nyquistMask)
         % Detect near-zero phase points
        if ~isempty(find(abs(correctedWave(s,:)) < zeroPhaseThreshold))
                candidates = find(abs(correctedWave(s,:)) < zeroPhaseThreshold);
                if numel(candidates) > 1
                    % Generation of variants
                    variants = cell(numel(candidates), 1);
                    scores = zeros(numel(candidates), 1);

                    for c = 1:numel(candidates)
                        variant = correctedWave(s,:);
                        anchor = candidates(c);
                        
                        % Left Correction 
                        for a = anchor-1:-1:1
                            delta = variant(a+1) - variant(a);
                            if abs(delta) > phaseJumpThreshold
                                variant(a) = variant(a) + round(delta/(2*pi))*2*pi;
                            end
                        end
                        
                        % Right Correction 
                        for a = anchor+1:nAngles
                            delta = variant(a) - variant(a-1);
                            if abs(delta) > phaseJumpThreshold
                                variant(a) = variant(a) - round(delta/(2*pi))*2*pi;
                            end
                        end
                        
                        % Calculate score as RMSE from trusted profile (NEW)
                        scores(c) = sqrt(mean((variant - trustedProfile).^2));
                        variants{c} = variant;
                    end

                    % Select variant with smallest distance to trusted profile
                    [~, bestIdx] = min(scores);
                    correctedWave(s,:) = variants{bestIdx};

                else

                    [~, anchorIdx] = min(abs(correctedWave(s,:)));  % Finds global minimum index
        
                    % Correct leftward from anchor
                    for a = (anchorIdx-1):-1:1
                        delta = correctedWave(s,a) - correctedWave(s,a+1);
                        if delta > pi
                            correctedWave(s,a) = correctedWave(s,a) - 2*pi;
                        elseif delta < -pi
                            correctedWave(s,a) = correctedWave(s,a) + 2*pi;
                        end
                    end
        
                    % Correct rightward from anchor
                    for a = (anchorIdx+1):nAngles
                        delta = correctedWave(s,a) - correctedWave(s,a-1);
                        if delta > pi
                            correctedWave(s,a) = correctedWave(s,a) - 2*pi;
                        elseif delta < -pi
                            correctedWave(s,a) = correctedWave(s,a) + 2*pi;
                        end
                    end

                end
        else   
            %%% cissors method %%%
            % First pass: left → right
            for a = 2:nAngles
                delta = correctedWave(s,a) - correctedWave(s,a-1);
                if abs(delta) > phaseJumpThreshold
                    if delta > pi
                        correctedWave(s,a) = correctedWave(s,a) - 2*pi;
                    elseif delta < -pi
                        correctedWave(s,a) = correctedWave(s,a) + 2*pi;
                    end
                end
            end
    
            % second pass: right → left (compensation)
            for a = (nAngles-1):-1:1
                delta = correctedWave(s,a) - correctedWave(s,a+1);
                if abs(delta) > phaseJumpThreshold
                    if delta > pi
                        correctedWave(s,a) = correctedWave(s,a) - 2*pi;
                    elseif delta < -pi
                        correctedWave(s,a) = correctedWave(s,a) + 2*pi;
                    end
                end
            end
    
            %%% Result merging %%%
            % Weighted Mean of the two passes 
            correctedWave(s,:) = 0.5*(correctedWave(s,:) + flip(correctedWave(s,:)));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 2: Constrained Trend Smoothing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trustedMask = (Scales > max(nyquistScales)) & (Scales <= upperCutoff);
    %%
    if any(trustedMask)
        % Calculate baseline trend with scale-dependent decay
        if matrixMode
            baseline = zeros(nScales, nAngles);
            % Retrieve the 'calibrationMatrix' from base workspace
            Matrix = evalin('base','calibrationMatrix'); 

            for s = 1:nScales
                thisScale = Scales(s);
        
                % Interpolate factor from the user-provided table
                factor = interp1( ...
                    Matrix (:,1), ...  % known scale values
                    Matrix (:,2), ...  % known factors
                    thisScale,           ...    % the scale we have
                    'linear', 'extrap' );
        
                baseline(s,:) = mean(correctedWave(trustedMask,:), 1, 'omitnan') * factor;
            end

        else
            % Use an exponential factor
            baseline = mean(correctedWave(trustedMask,:), 1, 'omitnan') ...
                       .* exp(-decaySharpness * (Scales'/max(Scales))) * decayFactor;
        end

    else
        baseline = zeros(nScales, 1);
    end

    % Apply constrained smoothing for large scales
    for s = find(Scales > upperCutoff)'
        correctedWave(s,:) = alpha*baseline(s,:) + (1-alpha)*correctedWave(s,:);
    end

    roiSpeedCell_final{iR} = correctedWave;
end
end

%==========================================================================

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

function createVideoFromFrames( ...
    fNames, fTimes, dataDir, varName, ...
    outputVideoFile, frameRate, ...
    applyShrink, shrinkfactor, invshrinkfactor, ...
    applyWindow, windowType, radius_factor, decay_rate, ...
    instrument, methodName, ...
    IR_threshold, IR_fillPercentile, ...
    VIS_lowerPercentile, VIS_upperPercentile, ...
    VIS_fillPercentile, clipMinHP, clipMaxHP, ...
    lpWidth20, lpWidth50, lpWidth100, ...
    Insolation_Correction, ...
    syntheticWaveMode, ...
    driftMode, ...
    time_resolution, ...
    drift_speed_m_s, driftAngleDeg, ...
    cphase, wavelength, direction, zamplitude, ...
    PBLdepth, dB_dzPBL, ...
    packet_center_x, packet_center_y, ...
    packet_width_x, packet_width_y, ...
    DXFactor)
%CREATEVIDEOFROMFRAMES Creates a video from a sequence of frames,
% applying the *same* modifications (synthetic waves, drift, preprocessing)
% that happen in the main processing loop, so the video matches final processed frames.
%
% Note: We assume:
%   - base_frame is the image from the FIRST file or some reference image
%     used for syntheticWaveMode/drift. If you used "base_frame = data from first frame"
%     in the main loop, pass that same matrix here.
%   - The time resolution between frames is time_resolution (seconds).
%   - Synthetic wave injection and drift logic are exactly as in your main loop.
%

numFrames = numel(fNames);
if numFrames == 0
    fprintf('No frames provided to create video.\n');
    return;
end

fprintf('Creating video: %s\n', outputVideoFile);

% --- Initialize Video Writer ---
try
    writerObj = VideoWriter(outputVideoFile, 'MPEG-4'); % Using MPEG-4 for compatibility
    writerObj.FrameRate = frameRate;
    open(writerObj);
catch ME
    fprintf('Error initializing VideoWriter: %s\n', ME.message);
    fprintf('Video creation failed.\n');
    return;
end

% For convenience, store original pixel size in km:
original_px_km = evalin('base','original_px_km');

% Loop over frames
for f_idx = 1:numFrames
    thisFileName = fNames{f_idx};
    thisFullPath = fullfile(dataDir, thisFileName);
    thisTime = fTimes(f_idx);
    fprintf('  Adding frame %d/%d: %s\n', f_idx, numFrames, thisFileName);

    % -- 1) Read raw data --
    try
        data = double(ncread(thisFullPath, varName));
    catch readME
        fprintf('  Warning: Could not read frame %s. Skipping. Error: %s\n', thisFileName, readME.message);
        continue; % Skip to next
    end

    % -- 2) If we do synthetic waves or drift, replicate main-loop logic --
    %    If f_idx == 1, we store data as base_frame in the main code
    %    but here we already have it passed in as "base_frame".
    %    So for f_idx>=2, override "data" with "base_frame" if synthetic is on.
    if f_idx == 1
        % We already have base_frame externally. So for consistency:
        if syntheticWaveMode || driftMode
            % We typically keep the original data for f_idx=1,
            base_frame = data;
            data = base_frame;
        end
    elseif f_idx ~= 1
        % If syntheticWaveMode or driftMode => we re-apply them
        if syntheticWaveMode || driftMode
            data = base_frame;  % exactly as in your main loop
        end
    end

    % => Synthetic wave injection
    if syntheticWaveMode
        % The same injection code from the main loop:
        [rowsF, colsF] = size(base_frame);
        [X, Y] = meshgrid(1:colsF, 1:rowsF);

        DX = 1000 * original_px_km * DXFactor;  % pixel size in meters
        Xm = (X - mean(X(:))) * DX;
        Ym = (Y - mean(Y(:))) * DX * -1;  % invert if y down

        % --- meteo → math angle conversion (0° = Est, trigonometrical direction) ---
        theta = deg2rad(90 - direction);    % direction given in meteorological convention
        
        % --- wave vector ---
        k = (2*pi / wavelength) * cos(theta);   % kx
        l = (2*pi / wavelength) * sin(theta);   % ky  (kept for phase)
        omega = cphase * (2 * pi / wavelength);

        t = (f_idx - 1) * time_resolution;  % e.g. in seconds
        phase = k * Xm + l * Ym - omega * t;
        
        % k     = (2 * pi / wavelength) * cosd(direction);
        % l     = (2 * pi / wavelength) * sind(direction);
        % omega = cphase * (2 * pi / wavelength);
        % 
        % t = (f_idx - 1) * time_resolution;  % e.g. in seconds
        % phase = k * Xm + l * Ym - omega * t;

        dz = zamplitude * sin(phase);

        % Envelope
        Ampwindow = exp( -(((Xm - packet_center_x) / packet_width_x).^2 ...
                        + ((Ym - packet_center_y) / packet_width_y).^2) );
        dz = dz .* Ampwindow;

        % --- horizontal displacement associayted to w' ---
        dxy = (zamplitude / PBLdepth) * wavelength .* ...
               sin(phase - pi/2) ./ DX;          % scalar amplitude in px
        dx  = dxy .* cos(theta).* Ampwindow;                 % composante O-E  (columns)
        dy  = dxy .* sin(theta).* Ampwindow;                 % composante S-N  (lines)


        % Horizontal displacement
        % dxy = (zamplitude / PBLdepth) * wavelength * sin(phase - pi/2) / DX;
        % dx = dxy .* cosd(direction) .* Ampwindow;
        % dy = dxy .* sind(direction) .* Ampwindow;

        XI = X - dx;
        YI = Y - dy;

        warped_img = interp2(X, Y, base_frame, XI, YI, 'linear', 0);
        modulated_img = warped_img .* (1 + dz / PBLdepth * dB_dzPBL);

        data = modulated_img;  % final synthetic wave injection
    end

    % => Drift shift
    if driftMode && f_idx > 1
        driftDistance_m = drift_speed_m_s * time_resolution * (f_idx - 1);
        driftDistance_km = driftDistance_m / 1000;
        pxShift = driftDistance_km / (original_px_km);  % approximate integer shift
        pxShift = round(pxShift);

        shift_dx =  pxShift * cosd(driftAngleDeg);
        shift_dy = -pxShift * sind(driftAngleDeg);
        shift_dx = round(shift_dx);
        shift_dy = round(shift_dy);

        data = circshift(data, [shift_dy, shift_dx]);
    end

    % -- 3) Preprocess frame (identical to main script) --
    data_pre = preprocessFrame(data, instrument, methodName, ...
                thisFullPath, thisTime, ...
                IR_threshold, IR_fillPercentile, ...
                VIS_lowerPercentile, VIS_upperPercentile, ...
                VIS_fillPercentile, clipMinHP, clipMaxHP, ...
                lpWidth20, lpWidth50, lpWidth100, ...
                Insolation_Correction);

    % -- 4) Optional shrinking/windowing for the video --
    if applyShrink && shrinkfactor ~= 1
        data_pre = imresize(data_pre, invshrinkfactor);
    end

    if applyWindow
        switch lower(windowType)
            case 'radial'
                data_pre = applyRadialWindow(data_pre, radius_factor, decay_rate);
            case 'rectangular'
                data_pre = applyRectangularWindow(data_pre, radius_factor, decay_rate);
        end
    end

    % Optionally invert IR
    if strcmpi(instrument,'IR')
        data_pre = -data_pre;
    end

    % -- 5) Convert final preprocessed image to a video frame --
    minVal = min(data_pre(:));
    maxVal = max(data_pre(:));
    if maxVal > minVal
        frame_norm = (data_pre - minVal) / (maxVal - minVal);
    else
        frame_norm = zeros(size(data_pre), 'like', data_pre); % handle constant image
    end

    frame_uint8 = uint8(frame_norm * 255);
    frame_rgb   = cat(3, frame_uint8, frame_uint8, frame_uint8); % grayscale RGB

    % -- 6) Write the frame to video --
    try
        writeVideo(writerObj, frame_rgb);
    catch writeME
        fprintf('  Warning: Could not write frame %d. Skipping. Error: %s\n', f_idx, writeME.message);
        % Optionally close the writer if repeated errors occur
    end
end

% --- Finalize Video ---
try
    close(writerObj);
    fprintf('Video creation complete: %s\n', outputVideoFile);
catch closeME
    fprintf('Error closing VideoWriter: %s\n', closeME.message);
end
end

function newName = makeParallelName(origFullPath, tag)
%   /path/IR_20231011_1400.nc  + '.spectralpeaks'  -->
%   /path/IR_20231011_1400.spectralpeaks.nc
[pathstr,base,~] = fileparts(origFullPath);
newName = fullfile(pathstr,[base tag '.nc']);
end

%------------------------------------------------------------------
function prof = azimProfile(img, cx, cy, radBins)
% AZIMPROFILE Returns the average radial profile ⟨img⟩(r)
%
% prof = azimProfile(img, cx, cy, radBins)
%
% img: 2-D image (NaN allowed)
% cx, cy: center (columns, rows) in pixels
% radBins: ring bounds [nBins+1 x 1]
%
% prof: vector 1×nBins containing azimuth average
%
% Notes
% -----
% - Ignores NaN with mean(...,'omitnan')
% - Handles cases where no pixel is present in the ring

    [Ny, Nx] = size(img);
    [X, Y]   = meshgrid(1:Nx, 1:Ny);
    R        = hypot(X - cx, Y - cy);
    

    nBins = numel(radBins)-1;
    prof  = NaN(1, nBins);

    for k = 1:nBins
        mask = (R >= radBins(k)) & (R < radBins(k+1));
        if any(mask(:))
            prof(k) = mean(img(mask), 'omitnan');
        end
    end
end
