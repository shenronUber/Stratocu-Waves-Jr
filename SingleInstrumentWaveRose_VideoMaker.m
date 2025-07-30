%% PASS-2  ·  PER-FRAME FILTERED ARRAYS  +  VIDEO WITH DYNAMIC CONTOURS
%     ––– relies on peakMaskROI saved in §6 ––– 
%
%  ⮑  For every raw frame it will …
%      1.  repeat the *same* preprocessing (synthetic wave, drift, high-pass …)
%      2.  run the CWT again
%      3.  for each ROI apply its own  (scale,angle)  mask  →  F_roi(x,y)
%      4.  write one NetCDF   IR_YYYYMMDD_HHMM.filtered.nc   (Ny × Nx × nROI)
%      5.  draw crest (red) / trough (blue) perimeters on the movie frame
%      6.  append the frame to  <instrument>_annot_video_<dateTag>.mp4
%
%  NOTE  – this block re-uses many variables already in memory
%          (fNames, fTimes, Scales, Angles, squares, …)

%% 1) VARIABLES SETUP
%----------- DATA-SOURCE MODE --------------------------------------------
useCustomFolder        = false;   % false ➜ normal GOES workflow
customFolderPath       = 'C:\Users\admin\Box\GWaves_Synthetic_G16ncfiles\closedcell';
customFilePattern      = 'closedcell_IR_2waves_halfhour%d.nc';   % sprintf() pattern
customFileIndices      = 0:9;          % which synthetic files to load
customStartDate        = datetime(2024,1,1,0,0,0);  % anchor timestamp of frame #0
custom_t_seconds       = 1800;         % spacing between successive frames


if useCustomFolder
    startDate = customStartDate;
    endDate = customStartDate + seconds(custom_t_seconds)*length(customFileIndices);
else
%----------- DATE/TIME SETTINGS ---------------------------
%startDate = datetime(2023, 10, 12, 1, 0, 0); % start processing time
%endDate   = datetime(2023, 10, 12, 3, 30, 0); % end processing time
startDate = datetime(2023, 10, 12, 14, 0, 0); % start processing time
%endDate = datetime(2023, 10, 11, 15, 0, 0); % start processing time
endDate = datetime(2023, 10, 12, 20, 30, 0); % end processing time
end

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
IR_threshold = 0; %277;         % IR threshold: values below are set to NaN
IR_fillPercentile = 50;     % Fill IR masked pixels with this percentile
useCumulativeMask   = true; % If true we build one static mask
cumulativeIRmask    = [];   % Will hold the union of clouds
True_Color_IR = true;       % Will display the IR video in real colors ( harder to see waves but real visuals )

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

%----------- SYNTHETIC DATA SETTINGS ----------------------
syntheticWaveMode = false;   % If true, superimpose a synthetic wave on a fixed base image
driftMode         = false;   % If true, apply a drift shift each frame using circshift

% Drift parameters (in m/s) rather than pixels/frame:
drift_speed_m_s = 15;       % e.g. 10 m/s
driftAngleDeg   = 45+90;       % e.g. 45 degrees (0 = right, 90 = up) // the image is inverted !!


% Parameters for the synthetic wave:
cphase = 15;                   % Phase speed (m/s)
wavelength = 2*150e3;            % Wavelength in meters
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

singleOutDir = 'C:\Users\admin\Box\GWaves_2023_10_11-14_SEPAC\videos';
if ~exist(singleOutDir, 'dir')
    mkdir(singleOutDir);
end

%----------- VIDEO OUTPUT SETTINGS -----------------------
createOutputVideo = true;    % Set to true to generate a video of processed frames
videoFrameRate = 5;          % Frames per second for the output video
videoApplyShrink = false;     % Use true to apply shrinkfactor to video frames
videoApplyWindow = true;     % Use true to apply windowing to video frames

overlayWaveContoursOnVideo = false;   % draw crest/trough contours ?
lineW       = 2;      % cosmetic

%----------- PEAK-DETECTION/Brightness decay SETTINGS -----------------------------------
contourOption       = 'percentile';           % 'percentile' | '3sigma'
contourArray        = [50 60 70];             % if 'percentile'


%% 2) RETRIEVE FILE LIST --------------------------------------------------
if useCustomFolder
    % -- synthetic experiment --------------------------------------------
    dataDir = customFolderPath;                     % skip <rootSepacDir>\IR\Data
    fNames  = arrayfun(@(k) sprintf(customFilePattern,k), ...
                       customFileIndices, 'uni',0);
    
    % forge deterministic timestamps so dateTag & filenames stay consistent
    fTimes  = customStartDate + seconds(custom_t_seconds)*(0:numel(fNames)-1).';
    startDate = fTimes(1);              % keep the rest of the script happy
    %endDate   = fTimes(end);            %   (used only for pretty strings)

    % detect which variable to read only once
    info      = ncinfo(fullfile(dataDir,fNames{1}));
    varList   = {info.Variables.Name};
    if  any(strcmp(varList,'CMI')),  varName = 'CMI';
    elseif any(strcmp(varList,'Rad')), varName = 'Rad';
    else  varName = varList{1}; 
    end
else
    % -- normal GOES hierarchy -------------------------------------------
    dataDir  = fullfile(rootSepacDir,upper(instrument),'Data');
    if ~exist(dataDir,'dir')
        error('Data directory for %s not found: %s',instrument,dataDir);
    end
    [fNames, fTimes, varName] = getDateRangeFiles(dataDir,startDate,endDate);
end

numFrames = numel(fTimes);
if numFrames==0
    error('No frames found for %s in the selected mode.',instrument);
end

fprintf('Found %d frames for instrument %s.\n',numFrames,instrument);


if useCumulativeMask
    fprintf('\n--- Pass-0 : building cumulative high-cloud mask (drift-aware) ---\n');
    
    cumulativeIRmask = false;          % au format après shrink & window
    for i = 1:numFrames
        % 1) raw mask from the original (un-drifted) frame
        data = double( ncread( fullfile(dataDir,fNames{i}), varName ) );
        
        % 2) apply transforms to the image
    
        % ---- On the first frame, store base_frame & possibly precompute synthetic-wave grids ----
        if i == 1
            base_frame = data;
    
            [rowsF, colsF] = size(base_frame);
            if syntheticWaveMode                
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
            t = (i - 1) * time_resolution;  % e.g. in seconds
    
            
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
        if driftMode && i > 1
            % Determine how many meters the image should shift in the time between frames:
            driftDistance_m = drift_speed_m_s * time_resolution* (i - 1);  % e.g. 10 m/s * 1800 s = 18000 m
            
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
    
        thisMask = data < IR_threshold;        % high clouds in *this* frame
    
        % 3) logical union
        if i == 1
            cumulativeIRmask = thisMask;
        else
            cumulativeIRmask = cumulativeIRmask | thisMask;
        end
    end
    fprintf('☑  cumulative mask created (%d×%d)\n', size(cumulativeIRmask));
end

%% Build ROI squares.
temp=round(rowsF/shrinkfactor);
rowsF = round(colsF/shrinkfactor);
colsF = temp;

% fprintf('Dimensions utilisées pour calculer la grille des ROI (après transposition et réduction):\n');
% fprintf('rowsF = %d, colsF = %d\n', rowsF, colsF);

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


%% Running the PASS-2

% 1) recreate the same dateTag you used when you saved it
dateTag = sprintf('%s_to_%s', ...
    datestr(startDate,'yyyymmddHHMM'), ...
    datestr(endDate,  'yyyymmddHHMM'));

% 2) build the full path
maskFile = fullfile(rootSepacDir,'PEAKMASK', ...
           sprintf('%s_peakmaskROI_%s.mat', instrument, dateTag));

% 3) load the variables you need
%    — peakMaskROI  : [NSCALES × NANGLES × nROI]
%    — Scales, Angles : your wavelet parameters
%    — squares      : the ROI grid (shrunken coordinates)
load(maskFile, 'peakMaskROI','Scales','Angles');

fprintf('\n──────────────── PASS-2  ────────────────\n');
dateTag   = sprintf('%s_to_%s', ...
              datestr(startDate,'yyyymmddHHMM'), ...
              datestr(endDate,  'yyyymmddHHMM'));

filteredDir = fullfile(rootSepacDir,'FILTERED');
if ~exist(filteredDir,'dir');  mkdir(filteredDir);  end

outVid = fullfile(singleOutDir, ...
          sprintf('%s_annot_video_%s.mp4', upper(instrument), dateTag));
vw = VideoWriter(outVid,'MPEG-4');
vw.FrameRate = videoFrameRate;
open(vw);

% handy constants
NS = NSCALES;      NA = NANGLES;    nROI=numSquares;  

thisFullPath = fullfile(dataDir, fNames{1});
data = double( ncread(thisFullPath, varName) );
% disp(size(data));

[NxFull, NyFull] = size(data);     

borderMaskFull = false(NyFull, NxFull);

for r = 1:numSquares
    % squares(r).x_range / y_range are in *shrunken* coordinates
    xShr = squares(r).x_range;
    yShr = squares(r).y_range;

    if shrinkfactor ~= 1
        xFull = (xShr(1)-1)*shrinkfactor + (1 : numel(xShr)*shrinkfactor);
        yFull = (yShr(1)-1)*shrinkfactor + (1 : numel(yShr)*shrinkfactor);
        % S'assurer de ne pas dépasser les bornes
        xFull = xFull( xFull <= NxFull );
        yFull = yFull( yFull <= NyFull );
    else
        xFull = xShr;
        yFull = yShr;
    end

    % draw the 1-pixel border
    borderMaskFull( yFull([1 end]), xFull            ) = true;  % top/bottom
    borderMaskFull( yFull,           xFull([1 end])  ) = true;  % left/right
end

% thicken with the same structuring element
se = strel('disk', lineW, 0);
borderMaskThick = imdilate(borderMaskFull, se);



for f_idx = 1:numFrames
    thisTime     = fTimes(f_idx);
    frameDateStr = datestr(thisTime,'yyyy_mm_dd_HHMMSS');
    fprintf('▶ Pass-2  frame  %d / %d  (%s)\n',f_idx,numFrames,frameDateStr);

    %% 7-A  read & preprocess exactly like §3  --------------------------
    thisFullPath = fullfile(dataDir, fNames{f_idx});
    data = double( ncread(thisFullPath, varName) );

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
                Insolation_Correction,cumulativeIRmask);
    
    % fprintf('Dimensions de l''image APRÈS preprocessFrame (contient la transposition):\n');
    % disp(size(data_pre));

    % ---- Resize or window if needed ----
    if shrinkfactor ~= 1
        data_pre = imresize(data_pre, invshrinkfactor); 
    end
    % fprintf('Dimensions FINALES de data_pre (après resize):\n');
    [Nx2, Ny2] = size(data_pre);

    % disp(size(data_pre));
    if doWindow
        switch lower(windowType)
            case 'radial'
                data_pre = applyRadialWindow(data_pre, radius_factor, decay_rate);
            case 'rectangular'
                data_pre = applyRectangularWindow(data_pre, radius_factor, decay_rate);
        end
    end

    %% 7-B  CWT & per-ROI filtered arrays ------------------------------
    waveStruct = cwtft2(data_pre,'wavelet','cauchy', ...
                        'scales',Scales,'angles',Angles);
    spec_full  = squeeze(waveStruct.cfs);
    for s=1:NS; spec_full(:,:,s,:) = spec_full(:,:,s,:) * ((pi/sqrt(2))/Scales(s)); end

    % allocate container  Ny × Nx × nROI
    Fstack = NaN(Nx2,Ny2,nROI,'single');

    for iROI = 1:nROI
        SA_mask = peakMaskROI(:,:,iROI);
        if ~any(SA_mask(:)); continue; end
        F = zeros(Nx2,Ny2,'single');
        [sIdx,aIdx] = find(SA_mask);
        for k = 1:numel(sIdx)
            F = F + spec_full(:,:, sIdx(k), aIdx(k));
        end
        Fstack(:,:,iROI) = F;
    end

    %% 7-C  crest / trough outlines for the video ----------------------
    Fstack = real(Fstack);

    if shrinkfactor ~= 1
    FstackFull = zeros([size(imresize(Fstack(:,:,1),shrinkfactor)) 30]);
    FstackFull = FstackFull(1:end-1,:,:);
    crestFrame  = false(size(imresize(Fstack(:,:,1),shrinkfactor)));
    crestFrame  = crestFrame(1:end-1,:,:);
    troughFrame = false(size(imresize(Fstack(:,:,1),shrinkfactor)));
    troughFrame = troughFrame(1:end-1,:,:);
    else
    FstackFull = zeros(size(Fstack));
    crestFrame  = false(Nx2,Ny2);
    troughFrame = false(Nx2,Ny2);
    end

    for iROI = 1:nROI
        wav_real = Fstack(:,:,iROI);
        
        % Restrict the display to the current ROI only ----
        xR       = squares(iROI).x_range;
        yR       = squares(iROI).y_range;
        mask = false(size(wav_real));
        % Mettre ce print juste avant la ligne qui plante
        % fprintf('iROI=%d: Tentative d''indexer la matrice mask (taille %dx%d) avec yR max=%d et xR max=%d\n', ...
        % iROI, size(mask,1), size(mask,2), max(yR), max(xR));

        mask(yR, xR) = true;

        % fprintf('Taille du masque APRÈS indexation: %dx%d\n', size(mask,1), size(mask,2));

        wav_real(~mask) = NaN;

        if all(isnan(wav_real(:))) || std(wav_real(:),'omitnan') == 0   % <-- skip ROI 
            continue
        end

        %=====  up-sample to full resolution  ======================
        if shrinkfactor ~= 1
            coeff_full  = imresize(wav_real,   shrinkfactor, 'bilinear');  % complex
            mask_full   = imresize(mask,    shrinkfactor, 'nearest')>0;
            coeff_full = coeff_full(1:end-1,:);
            mask_full = mask_full(1:end-1,:);
        else
            coeff_full  = wav_real;   % already full-res
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

        crest  = bwperim(crest);
        trough = bwperim(trough);
        
        % Merges with the global masks
        crestFrame  = crestFrame  | crest;
        troughFrame = troughFrame | trough;
             
    end
    %% 7-D  write per-frame NetCDF  (Nx × Ny × nROI) -------------------
    ncOut = fullfile(filteredDir, ...
            sprintf('%s_%s.filtered.nc', instrument,frameDateStr));
    if exist(ncOut,'file'); delete(ncOut); end

    if shrinkfactor ~= 1
    nccreate(ncOut,'FILTERED', 'Datatype','single', ...
                     'Dimensions',{'y',(Nx2*shrinkfactor)-1,'x',Ny2*shrinkfactor,'roi',nROI});
    else
    nccreate(ncOut,'FILTERED', 'Datatype','single', ...
                     'Dimensions',{'y',Nx2,'x',Ny2,'roi',nROI});     
    end

    ncwrite (ncOut,'FILTERED',FstackFull);

    %% 7-E  build the RGB video frame ----------------------------------
    % --- assure que data_pre et crestFrame sont de même taille -------------
    if shrinkfactor ~= 1    % crestFrame est Ny×Nx-1
        bg = imresize(data_pre, shrinkfactor, 'bilinear');
        bg = bg(1:size(crestFrame,1), :);          % retire éventuelle ligne fantôme
    else
        bg = data_pre;
    end
    data_pre = bg;

    minVal = min(data_pre(:));
    maxVal = max(data_pre(:));
    if maxVal > minVal
        frame_norm = (data_pre - minVal) / (maxVal - minVal);
    else
        frame_norm = zeros(size(data_pre), 'like', data_pre); % handle constant image
    end
    
    frame_uint8 = uint8(frame_norm * 255);
    frame_rgb   = cat(3, frame_uint8, frame_uint8, frame_uint8); % grayscale RGB

    se    = strel('disk', lineW, 0);  % disk ≃ circle
    
    % Dilate contours to desired thickness …
    crestDil  = imdilate(crestFrame , se);
    troughDil = imdilate(troughFrame, se);
    
    % --- candidate pixels that lie on a border ------------------------
    candCrest  = crestDil  & borderMaskThick;
    candTrough = troughDil & borderMaskThick;

    crestDil(candCrest)=0;
    troughDil(candTrough)=0;

    crestFrame = crestDil;
    troughFrame = troughDil;

    % paint
    R = frame_rgb(:,:,1); G = R; B = R;
    R(crestFrame)  = 255;  G(crestFrame)  = 0;  B(crestFrame)  = 0;
    R(troughFrame) = 0;    G(troughFrame) = 0;  B(troughFrame) = 255;
    frame_rgb = cat(3,R,G,B);     % l’intérieur reste gris

    writeVideo(vw, frame_rgb);

end % loop frames

close(vw);
fprintf('☑  Annotated video written → %s\n', outVid);
fprintf('☑  Per-frame filtered NetCDFs saved in  %s\n', filteredDir);



%% HELPER FUNCTIONS

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
    lpWidth20, lpWidth50, lpWidth100,Insolation_Correction,staticMask)
% PREPROCESSFRAME
%  Applies data-type-specific thresholds and then calls the chosen
%  method-based process. Fills or clamps as needed.

    switch upper(dataType)
        case 'IR'
            % (1) IR-specific threshold or masking
            if ~evalin( 'base', 'exist(''useCumulativeMask'',''var'') == 1' ) || isempty(staticMask)
                data(data < IR_threshold) = NaN; 
                nan_mask = isnan(data);     % frame-by-frame (old behaviour)
            elseif evalin( 'base', 'exist(''useCumulativeMask'',''var'') == 1' )
                nan_mask = staticMask;              % ⟵ NEW : same mask for everyone
            end
            data(nan_mask) = NaN;

                            
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
    instrument = evalin('base', 'instrument');
    True_Color_IR = evalin('base', 'True_Color_IR');

    if strcmpi(instrument, 'IR') && True_Color_IR
        % For IR, low data values (cold) are bright clouds.
        % Normal standardization would set these values to negative Z-scores.
        % By inverting the data before calculating the mean and standard deviation,
        % cold clouds will have positive Z-scores.
        data = -data;
        fprintf('Inverting IR data before standardization for physical consistency.\n');
    end
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
