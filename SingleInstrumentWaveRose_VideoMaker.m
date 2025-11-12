%% PASS-2  ·  TIERED FILTERING WITH MORPHOLOGICAL RECONSTRUCTION
%     ––– FINAL VERSION 3.0 –––
%
%  ⮑  For every raw frame PAIR it will …
%      1.  Preprocess two consecutive frames.
%      2.  Run the CWT on both frames to get the cross-spectrum.
%      3.  De-advect to find intrinsic phase/speed.
%      4.  Apply tiered filtering to find "MAYBE" and "CERTAIN" wave packets.
%      5.  Apply a 3-step morphological reconstruction to clean and connect crests/troughs.
%      6.  Write NetCDF/MAT files with tiered classification and ellipse data.
%      7.  Generate TWO annotated videos with reconstructed crest/trough contours.

%% 1) VARIABLES SETUP
%----------- DATA-SOURCE MODE --------------------------------------------
useCustomFolder        = true;   % false ➜ normal GOES workflow
customFolderPath       = 'C:\Users\admin\Box\GWaves_Synthetic_G16ncfiles\closedcellTest';
customFilePattern      = 'closedcell_IR_2waves_halfhour%d.nc';   % sprintf() pattern
customFileIndices      = 0:9;          % which synthetic files to load
customStartDate        = datetime(2104,1,1,0,0,0);  % anchor timestamp of frame #0
custom_t_seconds       = 1800;         % spacing between successive frames


if useCustomFolder
    startDate = customStartDate;
    endDate = customStartDate + seconds(custom_t_seconds)*length(customFileIndices);
else
%----------- DATE/TIME SETTINGS ---------------------------
startDate = datetime(2023, 10, 12, 14, 0, 0); % start processing time
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

pixel_size_km = original_px_km * shrinkfactor;

%----------- WAVELET PARAMETERS -----------------------------
Angles = 0 : pi/(7*2) : pi;           % Wavelet angles (in radians)
Scales = [2, 4, 8, 16, 32, 64];       % Wavelet scales in pixel units
NANGLES = numel(Angles);                  % Number of angles
NSCALES = numel(Scales);                  % Number of scales

%----------- WINDOWING SETTINGS -----------------------------
doWindow = true;              % Flag to apply windowing
windowType = 'rectangular';   % 'radial' or 'rectangular'
radius_factor = 0.6;          % Parameter for window function (if used)
decay_rate = 10;              % Controls steepness of window edge

%----------- PREPROCESSING THRESHOLDS -----------------------
if useCustomFolder
    IR_threshold = 0;    
else
    IR_threshold = 274; %274;
end
IR_fillPercentile = 50;
useCumulativeMask   = true;
cumulativeIRmask    = [];
True_Color_IR = true;
VIS_lowerPercentile = 10;
VIS_upperPercentile = 99;
VIS_fillPercentile = 50;
Insolation_Correction = true;

%----------- HIGH-PASS FILTER SETTINGS ----------------------
clipMinHP = -3; clipMaxHP = 3;
lowPassFilterWidth_50 = 50;

%----------- PREPROCESSING METHOD SELECTION -----------------
switch upper(instrument)
    case 'IR'; methodName = 'highpass_50_sqrt';
    case 'VIS'; methodName = 'none';
    otherwise; error('Unknown instrument: %s', instrument);
end

%----------- SYNTHETIC DATA SETTINGS ----------------------
syntheticWaveMode = false; time_resolution = 1800;

%----------- VIDEO & CONTOUR SETTINGS -----------------------
singleOutDir = 'C:\Users\admin\Box\GWaves_2023_10_11-14_SEPAC\videos';
if ~exist(singleOutDir, 'dir'); mkdir(singleOutDir); end
lineW = 2; % contour line width
contour_percentile = 70; % Use this percentile of real-part amplitude for contours
min_crest_trough_size_px = 150; % Minimum area for a region to be considered a "core"

% --- Saliency + persistence (global params) ---
saliency_decay          = 0.85; %0.85 or 1;   % EMA decay (0.85–0.95)
saliency_threshold      = 0.4;    % Threshold before considered vanished
min_region_area_px      = 300;    % drop tiny specks
saliency_se_disk        = 2;      % small gap closing

% persistence accumulator (created on first iteration)
saliency_accum = [];

% --- NEW PARAMETER for Reconstruction ---
morph_closing_radius_px = 10; % Radius (in pixels) for bridging gaps. Increase to connect more distant fragments.

% per-frame annotation cache + science accumulator + records 
annotDir = fullfile(rootSepacDir, 'ANNOTATION_FRAMES');
if ~exist(annotDir,'dir'); mkdir(annotDir); end

save_annotation_as_fullres = true;   % keep true to cut at final resolution
annotation_dtype = 'single';         % 'single' (fast) or 'int16' with scaling

%  CWT cache + Phase-B parameters 
cwtCacheDir = fullfile(rootSepacDir, 'CWT_CACHE', upper(instrument));
if ~exist(cwtCacheDir,'dir'); mkdir(cwtCacheDir); end
saveCWTcache      = true;    % write per-frame CWTs during Phase A
reuseCWTcache     = true;    % skip recomputing if cache present
deleteCWTcacheEnd = true;   % set true to delete after video

finalZonesFile = fullfile(rootSepacDir, sprintf('%s_final_zones_%s.mat', ...
                 upper(instrument), datestr(startDate,'yyyymmddHHMM')));

% Phase-B selection thresholds
meanPowerPercentile = 85;   % per-ROI mean power rose threshold
RbarMin             = 0.35; % circular-mean resultant length threshold
coh2Min             = 0.20; % min coherence^2 on time-mean cross
vmin_mps            = 7.0;  % min |speed| for blob

%% 2) RETRIEVE FILE LIST --------------------------------------------------
if useCustomFolder
    dataDir = customFolderPath;
    fNames  = arrayfun(@(k) sprintf(customFilePattern,k), customFileIndices, 'uni',0);
    fTimes  = customStartDate + seconds(custom_t_seconds)*(0:numel(fNames)-1).';
    startDate = fTimes(1);
    info = ncinfo(fullfile(dataDir,fNames{1}));
    varList = {info.Variables.Name};
    if any(strcmp(varList,'CMI')); varName = 'CMI';
    elseif any(strcmp(varList,'Rad')); varName = 'Rad';
    else; varName = varList{1}; end
else
    dataDir  = fullfile(rootSepacDir,upper(instrument),'Data');
    if ~exist(dataDir,'dir'); error('Data directory for %s not found: %s',instrument,dataDir); end
    [fNames, fTimes, varName] = getDateRangeFiles(dataDir,startDate,endDate);
end
numFrames = numel(fTimes);
if numFrames==0; error('No frames found for %s in the selected mode.',instrument); end
fprintf('Found %d frames for instrument %s.\n',numFrames,instrument);

if useCumulativeMask
    fprintf('\n--- Pass-0 : building cumulative high-cloud mask (drift-aware) ---\n');
    cumulativeIRmask = false;
    for i = 1:numFrames
        data = double( ncread( fullfile(dataDir,fNames{i}), varName ) );
        thisMask = data < IR_threshold;
        if i == 1; cumulativeIRmask = thisMask;
        else; cumulativeIRmask = cumulativeIRmask | thisMask; end
    end
    fprintf('☑  cumulative mask created (%d×%d)\n', size(cumulativeIRmask));
end

%% 3) LOAD ADVECTION DATA FROM PASS-1
fprintf('Loading advection data from Pass-1...\n');
dateTag = sprintf('%s_to_%s', datestr(startDate,'yyyymmddHHMM'), datestr(endDate,  'yyyymmddHHMM'));
advectionDir = fullfile(rootSepacDir, 'ADVECTION_DATA');
advectionFile = fullfile(advectionDir, sprintf('%s_advection_vectors_%s.mat', instrument, dateTag));
if ~exist(advectionFile, 'file'); error('Advection data file not found: %s', advectionFile); end
adv_data = load(advectionFile);
fprintf('Advection data loaded successfully.\n');

%% 4) NEW PASS-2 PARAMETERS AND SETUP
% --- PASS-1 reuse (expressed on the PASS-2 scale axis) ---
beta = 0.8;         % Weight for smoothing toward baseline trend
decayFactor = 0.55;  % Amplitude decay applied to baseline smoothing
decaySharpness = 1.2;  % Controls how sharply the baseline trend decays with scale

nyquistScales_post = [2 4 8];   % OK already (original px)
upperCutoff_post   = 16;        % original px; 8 was post-shrink
scalesForAdvection = [4 8 16];  % DO NOT divide by shrinkfactor

matrixMode = true;
calibrationMatrix = [
     2,   0.50;
     4,   0.50;
     8,   0.50;
    16,   0.40;
    32,   0.335;
    64,   0.17;
   128,   0.083;
];


% ROI expansion for PASS-1-on-ROI
pass1ROI.expandFactor = 1.75;     % 1.5–2 works well
pass1ROI.minW         = 64;       % px, ensure at least a few wavelengths
pass1ROI.minH         = 64;

cohThreshold       = 0.2;

advection_interp_method = 'linear';

% --- TIERED FILTERING PARAMETERS (fine-tuned) ---
power_threshold_factor = 1.5;
filterTiers.MAYBE.intrinsic_speed_threshold_ms = 7.5;
filterTiers.MAYBE.min_packet_size_pixels       = 200;
filterTiers.MAYBE.eccentricity_threshold       = 0.70;
filterTiers.MAYBE.solidity_threshold           = 0.5;
filterTiers.CERTAIN.intrinsic_speed_threshold_ms = 10;
filterTiers.CERTAIN.min_packet_size_pixels       = 500;
filterTiers.CERTAIN.eccentricity_threshold       = 0.70;
filterTiers.CERTAIN.solidity_threshold           = 0.6;

wavePacketDir = fullfile(rootSepacDir, 'WAVE_PACKET_DATA_TIERED');
if ~exist(wavePacketDir, 'dir'); mkdir(wavePacketDir); end

firstFullPath = fullfile(dataDir, fNames{1});
lat_vec = ncread(firstFullPath, 'latitude');
lon_vec = ncread(firstFullPath, 'longitude');
[Ny_full, Nx_full] = size(ncread(firstFullPath, varName)');

%% 5) MAIN PROCESSING LOOP (FRAME-PAIR ANALYSIS)
fprintf('\n──────────────── NEW PASS-2: SPATIAL WAVE PACKET DETECTION ────────────────\n');

for f_idx = 1 : numFrames - 1
    
    thisTime = fTimes(f_idx);
    fprintf('▶ Processing frame pair starting at %s\n', datestr(thisTime));
    
    % --- 5A) Preprocess consecutive frames ---
    path1 = fullfile(dataDir, fNames{f_idx});
    data1_raw = double(ncread(path1, varName));
    data1_pre = preprocessFrame(data1_raw, instrument, methodName, path1, thisTime, IR_threshold, IR_fillPercentile, VIS_lowerPercentile, VIS_upperPercentile, VIS_fillPercentile, clipMinHP, clipMaxHP, lowPassFilterWidth_50, Insolation_Correction, cumulativeIRmask);
    if shrinkfactor ~= 1; data1_pre = imresize(data1_pre, invshrinkfactor); end
    if doWindow; data1_pre = applyRectangularWindow(data1_pre, radius_factor, decay_rate); end

    path2 = fullfile(dataDir, fNames{f_idx+1});
    data2_raw = double(ncread(path2, varName));
    data2_pre = preprocessFrame(data2_raw, instrument, methodName, path2, fTimes(f_idx+1), IR_threshold, IR_fillPercentile, VIS_lowerPercentile, VIS_upperPercentile, VIS_fillPercentile, clipMinHP, clipMaxHP, lowPassFilterWidth_50, Insolation_Correction, cumulativeIRmask);
    if shrinkfactor ~= 1; data2_pre = imresize(data2_pre, invshrinkfactor); end
    if doWindow; data2_pre = applyRectangularWindow(data2_pre, radius_factor, decay_rate); end

    
    % --- 5B) Compute or load per-frame CWTs (cache-aware) ---
    cwtPath1 = fullfile(cwtCacheDir, sprintf('cwt_f%06d.mat', f_idx));
    cwtPath2 = fullfile(cwtCacheDir, sprintf('cwt_f%06d.mat', f_idx+1));
       
    % frame f_idx
    if reuseCWTcache && exist(cwtPath1,'file')
        S1 = load(cwtPath1,'W','meta'); 
        W1 = S1.W; 
        meta1 = S1.meta;
    else
        W1 = computeCWTframe(data1_pre, Scales, Angles);
        meta1 = struct('frame_index',f_idx,'time',fTimes(f_idx), ...
            'lat_vec',lat_vec,'lon_vec',lon_vec,'shrinkfactor',shrinkfactor, ...
            'Scales',Scales,'Angles',Angles,'pixel_size_km',pixel_size_km, ...
            'time_resolution',time_resolution);
        if saveCWTcache
            W = W1; meta = meta1;  % ensure variables exist with saved names
            save(cwtPath1,'W','meta','-v7.3');
        end
    end
    
    % frame f_idx+1
    if reuseCWTcache && exist(cwtPath2,'file')
        S2 = load(cwtPath2,'W','meta'); 
        W2 = S2.W; 
        meta2 = S2.meta;
    else
        W2 = computeCWTframe(data2_pre, Scales, Angles);
        meta2 = struct('frame_index',f_idx+1,'time',fTimes(f_idx+1), ...
            'lat_vec',lat_vec,'lon_vec',lon_vec,'shrinkfactor',shrinkfactor, ...
            'Scales',Scales,'Angles',Angles,'pixel_size_km',pixel_size_km, ...
            'time_resolution',time_resolution);
        if saveCWTcache
            W = W2; meta = meta2;  %#ok<NASGU>
            save(cwtPath2,'W','meta','-v7.3');
        end
    end
   
    % cross & power for CURRENT pair → used only for saliency/persistence (ROI finding)
    cross_spec = W1 .* conj(W2);
    power_spec = abs(W1).^2;
    [Ny_sh, Nx_sh, ~, ~] = size(cross_spec);

    % --- 5C) Interpolate advection field ---
    fprintf('   Step 1: Interpolating advection field...\n');
    [LON_sh, LAT_sh] = meshgrid(linspace(lon_vec(1), lon_vec(end), Nx_sh), linspace(lat_vec(1), lat_vec(end), Ny_sh));
    F_U = scatteredInterpolant(adv_data.roi_center_lon, adv_data.roi_center_lat, adv_data.advection_U, advection_interp_method, 'none');
    F_V = scatteredInterpolant(adv_data.roi_center_lon, adv_data.roi_center_lat, adv_data.advection_V, advection_interp_method, 'none');
    U_advection_sh = F_U(double(LON_sh), double(LAT_sh));
    V_advection_sh = F_V(double(LON_sh), double(LAT_sh));
    
    % --- 5D) Global de-advection used only for coarse detection ---
    fprintf('   Step 2: De-advecting each pixel...\n');
    observed_phase = angle(cross_spec);
    advection_phase = zeros(size(observed_phase), 'single');
    for iS = 1:NSCALES;lambda_m = Scales(iS) * (pixel_size_km*1000) * (pi/sqrt(2)); 
        for iA=1:NANGLES; adv_speed_in_dir=U_advection_sh*cos(Angles(iA))+V_advection_sh*sin(Angles(iA)); 
            advection_phase(:,:,iS,iA)=((adv_speed_in_dir*time_resolution)/lambda_m)*2*pi; 
        end 
    end
    intrinsic_phase = angle(exp(1i * (observed_phase - advection_phase)));
      
    % --- 5E) FIRST PASS: Initial packet detection using noisy pixel-level speed ---
    fprintf('   Step 3: Finding initial packet locations (first pass)...\n');
    dominant_intrinsic_speed_initial = zeros(Ny_sh, Nx_sh, 'single');
    power_threshold = mean(power_spec(:)) + power_threshold_factor * std(power_spec(:));
    for r = 1:Ny_sh; for c = 1:Nx_sh; [max_power, max_idx] = max(squeeze(power_spec(r, c, :, :)), [], 'all'); if max_power > power_threshold; [s_idx, a_idx] = ind2sub([NSCALES, NANGLES], max_idx); lambda_m_dom = (Scales(s_idx) * pixel_size_km * 1000) * pi / sqrt(2); dominant_intrinsic_speed_initial(r, c) = (intrinsic_phase(r, c, s_idx, a_idx) * lambda_m_dom) / (2 * pi * time_resolution); end; end; end
    
    % 5F) Salient & persistent region detection on TOP1 speed (relative, not absolute)

    V  = dominant_intrinsic_speed_initial;     % signed m/s
    Vm = abs(V);
    
    % --- A) Local z-score saliency of |V| (relative contrast, not absolute)
    h = fspecial('gaussian', 21, 5);                    % local window ~21×21
    mu_loc = imfilter(Vm, h, 'replicate');
    sd_loc = sqrt( max(imfilter(Vm.^2, h, 'replicate') - mu_loc.^2, 1e-6) );
    sd_loc(sd_loc<0.1) = NaN;

    % Flatten Z into a vector and compute threshold at desired quartile
    q = 0.80;                                % choose quartile: 0.75 = 75th percentile
    thresh = quantile(sd_loc(:), q);              % threshold value
    sd_loc(sd_loc < thresh) = 0;                % or NaN if you prefer
    
    % --- D) Hysteresis thresholding (quantile based, relative to current frame)
    T_hi = prctile(sd_loc(:), 98);     % seeds = top ~2%
    T_lo = prctile(sd_loc(:), 90);     % candidates = top ~10%
    seeds = sd_loc >= T_hi;
    cands = sd_loc >= T_lo;
    
    % Region growing via morphological reconstruction
    mask_salient = imreconstruct(seeds, cands);
    
    % Small cleanup
    se_small = strel('disk', saliency_se_disk);
    mask_salient = imclose(bwareaopen(mask_salient, min_region_area_px), se_small);

    % --- E) Persistence across frames (pixel-wise EMA)
    if isempty(saliency_accum), saliency_accum = zeros(size(mask_salient), 'single'); end
    saliency_accum = saliency_decay * saliency_accum + single(mask_salient);
    
    persistent_now = saliency_accum >= 3;           % present ≥ ~3 frames
    persistent_now = bwareaopen(persistent_now, min_region_area_px);     % drop tiny

end % End of main loop

    
% This is the mask we will study further
persistent_regions_mask = persistent_now;

%%  FINAL ZONES SUMMARY (ellipses + aggregated speeds)
assert(~isempty(persistent_regions_mask), 'persistent_regions_mask : no zone detected.');
science_bins = bwareaopen(persistent_regions_mask, 1000);

[Lz, Nz] = bwlabel(science_bins);
zoneProps = regionprops(Lz, 'Centroid','MajorAxisLength','MinorAxisLength','Orientation','Area');
fprintf('Zones finales trouvées : %d\n', Nz);

save(finalZonesFile, 'Lz','Nz','zoneProps','fTimes','Scales','Angles', ...
     'pixel_size_km','time_resolution','shrinkfactor','lat_vec','lon_vec','-v7.3');
fprintf('✓ Final zones saved to %s\n', finalZonesFile);

%% ─────────────────────────────  PHASE B (ex-Pass-3)  ─────────────────────────────
fprintf('\n──────────────── PASS-3: ROI spectral averaging & fixed-cell annotation ───────────────\n');

% Load final zones (already in memory, but make robust if rerun)
if ~exist('Lz','var') || ~exist('zoneProps','var')
    Sfz = load(finalZonesFile);
    Lz = Sfz.Lz; Nz = Sfz.Nz; zoneProps = Sfz.zoneProps;
    fTimes = Sfz.fTimes; Scales = Sfz.Scales; Angles = Sfz.Angles;
    pixel_size_km = Sfz.pixel_size_km; time_resolution = Sfz.time_resolution;
end
NS = numel(Scales); NA = numel(Angles);

% ROI masks on the CWT grid (same dims as W)
[NyZ, NxZ] = size(Lz);
roiMask = false(NyZ, NxZ, Nz);
for z=1:Nz
    roiMask(:,:,z) = ellipseMask(zoneProps(z), NyZ, NxZ);
    roiMask(:,:,z) = imerode(roiMask(:,:,z), strel('disk',1));  % shrink slightly
end

% --- DEBUG / DIAGNOSTIC CONTROLS FOR PASS-3 (running averages) ---
dbgRoiRose_enable   = true;         % turn ON/OFF running per-ROI rose figures
dbgRoiRose_zlist    = 1:min(3,Nz); %1:min(3,Nz);  % which zones to visualize (IDs)
dbgRoiRose_maxFrames= 13;            % show figures for first N frame pairs
dbgRoiRose_outDir   = fullfile(wavePacketDir,'ROI_ROSE_DEBUG');
if dbgRoiRose_enable && ~exist(dbgRoiRose_outDir,'dir'); mkdir(dbgRoiRose_outDir); end

% --- Build the summary figure on the first preprocessed frame
first_frame_path = fullfile(dataDir, fNames{1});
first_frame_raw  = double(ncread(first_frame_path, varName));
first_frame_pre  = preprocessFrame(first_frame_raw, instrument, methodName, ...
                                   first_frame_path, fTimes(1), IR_threshold, IR_fillPercentile, ...
                                   VIS_lowerPercentile, VIS_upperPercentile, VIS_fillPercentile, ...
                                   clipMinHP, clipMaxHP, lowPassFilterWidth_50, Insolation_Correction, ...
                                   cumulativeIRmask);
if shrinkfactor ~= 1; first_frame_pre = imresize(first_frame_pre, invshrinkfactor); end
if doWindow
    first_frame_pre = applyRectangularWindow(first_frame_pre, radius_factor, decay_rate);
end

% === PASS-3 accumulators (raw) ===
sumP1  = cell(Nz,1);   % ∑ mean power at t
sumP2  = cell(Nz,1);   % ∑ mean power at t+Δt
sumC   = cell(Nz,1);   % ∑ mean cross-spectrum (complex)
countF = zeros(Nz,1);  % number of contributing frame-pairs (scalar/ROI)

% Advection averaging (weighted) – used once at the end
sumU = zeros(Nz,1);    % ∑ (Ū_roi * weight)
sumV = zeros(Nz,1);    % ∑ (V̄_roi * weight)
sumW = zeros(Nz,1);    % ∑ weight

% Initialize shapes on first CWT
tmp = load(fullfile(cwtCacheDir, sprintf('cwt_f%06d.mat', 1)),'W');
[~,~,NS,NA] = size(tmp.W);

% Initialize accumulators with correct shapes
for z = 1:Nz
    sumP1{z} = zeros(NS,NA,'single');
    sumP2{z} = zeros(NS,NA,'single');
    sumC{z}  = zeros(NS,NA,'single');
end

% Iterate frame pairs again, reusing cached CWTs and advection de-advect
numFrames = numel(fTimes);
for f_idx = 1:numFrames-1
    % Load cached CWTs
    S1 = load(fullfile(cwtCacheDir, sprintf('cwt_f%06d.mat', f_idx)),   'W','meta');
    S2 = load(fullfile(cwtCacheDir, sprintf('cwt_f%06d.mat', f_idx+1)), 'W','meta');
    W1 = S1.W; W2 = S2.W;

    [Ny_sh, Nx_sh, NS, NA] = size(W1);

    cross_spec = W1 .* conj(W2);
    power1 = abs(W1).^2; 
    power2 = abs(W2).^2;

    % Interpolate advection to CWT grid
    [LON_sh, LAT_sh] = meshgrid(linspace(lon_vec(1), lon_vec(end), size(W1,2)), ...
                                linspace(lat_vec(1), lat_vec(end), size(W1,1)));
    F_U = scatteredInterpolant(adv_data.roi_center_lon, adv_data.roi_center_lat, ...
                               adv_data.advection_U, advection_interp_method, 'none');
    F_V = scatteredInterpolant(adv_data.roi_center_lon, adv_data.roi_center_lat, ...
                               adv_data.advection_V, advection_interp_method, 'none');
    U_advection_sh = F_U(double(LON_sh), double(LAT_sh));
    V_advection_sh = F_V(double(LON_sh), double(LAT_sh));

    for z=1:Nz
        M = roiMask(:,:,z); 
        np = nnz(M);
        if np==0, continue; end

        % Spatial means per ROI → [NS,NA]
        mP1 = squeeze(mean(reshape(power1(repmat(M,1,1,NS,NA)),[np,NS,NA]),1,'omitnan'));
        mP2 = squeeze(mean(reshape(power2(repmat(M,1,1,NS,NA)),[np,NS,NA]),1,'omitnan'));
        mC  = squeeze(mean(reshape(cross_spec(repmat(M,1,1,NS,NA)),[np,NS,NA]),1,'omitnan'));

        % Accumulate raw
        sumP1{z}  = sumP1{z} + mP1;
        sumP2{z}  = sumP2{z} + mP2;
        sumC{z}   = sumC{z}  + mC;
        countF(z) = countF(z) + 1;  % scalar

        % ROI-mean advection (for end-only de-advect) — use a robustness weight
        Ubar = mean(U_advection_sh(M),'omitnan'); 
        Vbar = mean(V_advection_sh(M),'omitnan');

        % Weight: total cross-power magnitude this frame in this ROI
        % (You could also use sum of mP1 or sqrt(mP1.*mP2); this one tracks coherence energy.)
        w_adv = nansum(abs(mC(:)));

        if isfinite(Ubar) && isfinite(Vbar) && isfinite(w_adv) && w_adv>0
            sumU(z) = sumU(z) + Ubar * w_adv;
            sumV(z) = sumV(z) + Vbar * w_adv;
            sumW(z) = sumW(z) + w_adv;
        end

        % ----------------- RUNNING VISUALIZATION (optional) -----------------
        if dbgRoiRose_enable && ismember(z, dbgRoiRose_zlist) && f_idx <= dbgRoiRose_maxFrames
            nFz = max(countF(z), 1);
        
            % running ensemble means (raw)
            meanP1_run = sumP1{z} ./ nFz;
            meanP2_run = sumP2{z} ./ nFz;
            meanC_run  = sumC{z}  ./ nFz;
        
            % --- provisional advection for display only ---
            if sumW(z) > 0
                Uhat = sumU(z) / sumW(z);
                Vhat = sumV(z) / sumW(z);
            else
                % fall back to this frame’s ROI-mean advection if no weights yet
                Uhat = Ubar; 
                Vhat = Vbar;
            end
        
            % build advection phase map (scale×angle)
            phi_adv_run = zeros(NS,NA,'single');
            for s_=1:NS
                lambda_m = Scales(s_) * (pixel_size_km*1000) * (pi/sqrt(2));
                for a_=1:NA
                    upar = Uhat*cos(Angles(a_)) + Vhat*sin(Angles(a_));
                    phi_adv_run(s_,a_) = ((upar * time_resolution)/lambda_m) * 2*pi;
                end
            end
        
            % observed → provisional intrinsic phase (running)
            phi_obs_run = angle(meanC_run);
            phi_int_run = angle(exp(1i*(phi_obs_run - phi_adv_run)));
        
            % (optional) apply the same Nyquist/scale cleaning just for the plot
            phi_int_cell_run = limitSpeedByScale({phi_int_run}, Scales, nyquistScales_post, upperCutoff_post, ...
                                                 beta, decayFactor, decaySharpness, true);
            phi_int_run = phi_int_cell_run{1};
        
            % diagnostics for running selection
            coh2_run = (abs(meanC_run).^2) ./ (meanP1_run .* meanP2_run);
        
            % phase concentration proxy at run time (since we already cleaned)
            Rbar_run = abs(exp(1i*phi_int_run));
        
            % convert phase → speed for masking and display
            v_run = zeros(NS,NA,'single');
            for s_=1:NS
                lambda_m = Scales(s_) * (pixel_size_km*1000) * (pi/sqrt(2));
                v_run(s_,:) = (phi_int_run(s_,:)/(2*pi)) .* (lambda_m/time_resolution);
            end
        
            % running mask (same criteria you use at the end)
            Pthr_run = prctile(meanP1_run(:), meanPowerPercentile);
            mask_run = (meanP1_run >= Pthr_run) & ...
                       (coh2_run >= coh2Min) & (abs(v_run) >= vmin_mps);
        
            % pick largest-power connected blob
            CCrun = bwconncomp(mask_run);
            msel_run = false(NS,NA);
            if CCrun.NumObjects>0
                best = 0; ib = 0;
                for kk=1:CCrun.NumObjects
                    idx = CCrun.PixelIdxList{kk};
                    sP  = sum(meanP1_run(idx),'omitnan');
                    if sP>best, best=sP; ib=kk; end
                end
                if ib>0, msel_run(CCrun.PixelIdxList{ib}) = true; end
            end
        
            % spatial panel uses your ROI mask (ellipse/intersection already baked in roiMask)
            roi_z = roiMask(:,:,z);
            bg_img = first_frame_pre;  % previously computed
        
            % plot three-panel running diagnostic
            % plotRunningRoiRose(z, f_idx, ...
            %     meanP1_run, meanP2_run, meanC_run, phi_int_run, msel_run, ...
            %     Scales, Angles, time_resolution, pixel_size_km, ...
            %     roi_z, bg_img, dbgRoiRose_outDir);
        end
    end
end
%% === END-ONLY de-advect + cleaning per ROI ===
% === Also includes Annotated video with fixed cells  ====
selectedCells  = cell(Nz,1);
speedMean_mps  = nan(Nz,1);
speedStd_mps   = nan(Nz,1);

for z = 1:Nz
    nF = max(countF(z), 1);

    % Ensemble means (raw)
    meanP1 = sumP1{z} ./ nF;
    meanP2 = sumP2{z} ./ nF;
    meanC  = sumC{z}  ./ nF;

    % Weighted mean advection (robust)
    if sumW(z) > 0
        Umean = sumU(z) / sumW(z);
        Vmean = sumV(z) / sumW(z);
    else
        Umean = 0; Vmean = 0;   % fallback
    end

    % Advection phase map (scale×angle)
    phi_adv = zeros(NS,NA,'single');
    for s = 1:NS
        lambda_m = Scales(s) * (pixel_size_km*1000) * (pi/sqrt(2));
        for a = 1:NA
            upar = Umean*cos(Angles(a)) + Vmean*sin(Angles(a));
            phi_adv(s,a) = ((upar * time_resolution)/lambda_m) * 2*pi;
        end
    end

    % Observed → intrinsic phase (single shot)
    phi_obs      = angle(meanC);
    phi_intr_raw = angle(exp(1i*(phi_obs - phi_adv)));

    % Nyquist/scale constraints (single shot)
    phi_cell     = limitSpeedByScale({phi_intr_raw}, Scales, nyquistScales_post, upperCutoff_post, ...
                                     beta, decayFactor, decaySharpness, true);
    mean_phi_int = phi_cell{1};  % radians

    % Phase → speed (m/s)
    mean_speed = zeros(NS,NA,'single');
    for s = 1:NS
        lambda_m = Scales(s) * (pixel_size_km*1000) * (pi/sqrt(2));
        mean_speed(s,:) = (mean_phi_int(s,:)/(2*pi)) .* (lambda_m/time_resolution);
    end

    % Coherence² for gating
    coh2 = (abs(meanC).^2) ./ (meanP1 .* meanP2);

    % Selection mask (power percentile + coherence + min speed)
    Pthr = prctile(meanP1(:), meanPowerPercentile);
    mask = (meanP1 >= Pthr) & (coh2 >= coh2Min) & (abs(mean_speed) >= vmin_mps);

    % Pick strongest connected blob by total mean power
    CC = bwconncomp(mask);
    msel = false(NS,NA);
    if CC.NumObjects > 0
        best = 0; ibest = 0;
        for k = 1:CC.NumObjects
            idx = CC.PixelIdxList{k};
            sP  = sum(meanP1(idx), 'omitnan');
            if sP > best, best = sP; ibest = k; end
        end
        msel(CC.PixelIdxList{ibest}) = true;
    end
    selectedCells{z} = msel;

    % Final ROI speed stats from the selected blob
    vvals = mean_speed(msel);
    speedMean_mps(z) = mean(vvals,'omitnan');
    speedStd_mps(z)  = std(vvals,'omitnan');
end

outVid = fullfile(singleOutDir, sprintf('%s_ZONE_ANNOT_PASS3_%s.mp4', ...
         upper(instrument), datestr(fTimes(1),'yyyymmddHHMM')));
vw = VideoWriter(outVid, 'MPEG-4'); vw.FrameRate = 10; open(vw);

for f_idx = 1:(numFrames-1)
    % Background image (reuse your preprocessFrame)
    path_bg = fullfile(dataDir, fNames{f_idx});
    raw_bg  = double(ncread(path_bg, varName));
    pre_bg  = preprocessFrame(raw_bg, instrument, methodName, ...
                              path_bg, fTimes(f_idx), IR_threshold, IR_fillPercentile, ...
                              VIS_lowerPercentile, VIS_upperPercentile, VIS_fillPercentile, ...
                              clipMinHP, clipMaxHP, lowPassFilterWidth_50, Insolation_Correction, ...
                              cumulativeIRmask);
    if doWindow, pre_bg = applyRectangularWindow(pre_bg, radius_factor, decay_rate); end
    bg = uint8(255 * mat2gray(pre_bg));  R=bg; G=bg; B=bg;

    % Load W for this frame
    S1 = load(fullfile(cwtCacheDir, sprintf('cwt_f%06d.mat', f_idx)),'W'); W1 = S1.W;
    [NyC, NxC, ~, ~] = size(W1);

    crest_all = false(NyC,NxC); trough_all = false(NyC,NxC);

    for z=1:Nz
        msel = selectedCells{z};
        if ~any(msel(:)), continue; end

        % reconstruct A from fixed cells
        A = zeros(NyC, NxC, 'single');
        [ss, aa] = ind2sub([NS,NA], find(msel));
        for k=1:numel(ss)
            A = A + real(W1(:,:,ss(k),aa(k)));
        end

        % ellipse mask at CWT resolution
        Ze = ellipseMask(zoneProps(z), NyC, NxC);
        vals = A(Ze);
        if isempty(vals) || all(~isfinite(vals)), continue; end
        t = prctile(abs(vals), contour_percentile);

        crest_all  = crest_all  | (Ze & (A >=  t));
        trough_all = trough_all | (Ze & (A <= -t));
    end

    se = strel('disk', lineW/2, 0);
    crest_perim  = imdilate(bwperim(crest_all),  se);
    trough_perim = imdilate(bwperim(trough_all), se);

    % Resize perimeters to bg if sizes differ
    if ~isequal(size(crest_perim), size(bg))
        crest_perim  = imresize(crest_perim,  size(bg), 'nearest');
        trough_perim = imresize(trough_perim, size(bg), 'nearest');
    end

    R(crest_perim)  = 255; G(crest_perim)  = 0;   B(crest_perim)  = 0;
    R(trough_perim) = 0;   G(trough_perim) = 0;   B(trough_perim) = 255;

    writeVideo(vw, cat(3,R,G,B));
end

close(vw);
fprintf('✓ PASS-3 video written → %s\n', outVid);

%% === FRAME-BY-FRAME v1-STYLE SPEEDS INSIDE PHASE B ===
% Outputs:
%   zone_frame_speeds{z}  -> [nFrames-1 x 1] mean speed per frame (m/s)
%   zone_frame_err{z}      -> [nFrames-1 x 1] std within-blob per frame (m/s)
%   speedMean_mps(z)       -> mean of per-frame means (m/s)
%   speedStd_mps(z)        -> std of per-frame means (random error, m/s)

zone_frame_speeds = cell(Nz,1);
zone_frame_err    = cell(Nz,1);

for z = 1:Nz
    zone_frame_speeds{z} = nan(numFrames-1,1);
    zone_frame_err{z}    = nan(numFrames-1,1);
end

for f_idx = 1:numFrames-1
    % load cached CWTs for this pair
    S1 = load(fullfile(cwtCacheDir, sprintf('cwt_f%06d.mat', f_idx)),   'W');
    S2 = load(fullfile(cwtCacheDir, sprintf('cwt_f%06d.mat', f_idx+1)), 'W');
    W1 = S1.W;  W2 = S2.W;                          % [Ny x Nx x NS x NA]
    cross_spec = W1 .* conj(W2);
    power1 = abs(W1).^2; 
    power2 = abs(W2).^2;
    [Ny_sh, Nx_sh, NS, NA] = size(W1);

    % interpolate advection to this CWT grid
    [LON_sh, LAT_sh] = meshgrid(linspace(lon_vec(1), lon_vec(end), Nx_sh), ...
                                linspace(lat_vec(1), lat_vec(end), Ny_sh));
    F_U = scatteredInterpolant(adv_data.roi_center_lon, adv_data.roi_center_lat, ...
                               adv_data.advection_U, advection_interp_method, 'none');
    F_V = scatteredInterpolant(adv_data.roi_center_lon, adv_data.roi_center_lat, ...
                               adv_data.advection_V, advection_interp_method, 'none');
    U_advection_sh = F_U(double(LON_sh), double(LAT_sh));
    V_advection_sh = F_V(double(LON_sh), double(LAT_sh));

    for z = 1:Nz
        M = roiMask(:,:,z); 
        np = nnz(M); 
        if np==0, continue; end

        % per-ROI spatial means (this frame pair) → [NS x NA]
        mP1 = squeeze(mean(reshape(power1(repmat(M,1,1,NS,NA)), [np,NS,NA]), 1, 'omitnan'));
        mP2 = squeeze(mean(reshape(power2(repmat(M,1,1,NS,NA)), [np,NS,NA]), 1, 'omitnan'));
        mC  = squeeze(mean(reshape(cross_spec(repmat(M,1,1,NS,NA)), [np,NS,NA]), 1, 'omitnan'));

        % frame-ROI mean advection
        Ubar = mean(U_advection_sh(M), 'omitnan');
        Vbar = mean(V_advection_sh(M), 'omitnan');
        if ~isfinite(Ubar), Ubar = 0; end
        if ~isfinite(Vbar), Vbar = 0; end

        % de-advect to intrinsic phase (v1 logic, per frame)
        phi_adv = zeros(NS,NA,'single');
        for s = 1:NS
            lambda_m = Scales(s) * (pixel_size_km*1000) * (pi/sqrt(2));
            for a = 1:NA
                upar = Ubar*cos(Angles(a)) + Vbar*sin(Angles(a));
                phi_adv(s,a) = ((upar * time_resolution)/lambda_m) * 2*pi;
            end
        end
        phi_obs  = angle(mC);
        phi_int  = angle(exp(1i*(phi_obs - phi_adv)));

        % optional: same Nyquist/scale cleaning used in v2
        phi_int = limitSpeedByScale({phi_int}, Scales, nyquistScales_post, upperCutoff_post, ...
                                    beta, decayFactor, decaySharpness, true);
        phi_int = phi_int{1};

        % phase → speed (m/s)
        v_map = zeros(NS,NA,'single');
        for s = 1:NS
            lambda_m = Scales(s) * (pixel_size_km*1000) * (pi/sqrt(2));
            v_map(s,:) = (phi_int(s,:)/(2*pi)) .* (lambda_m/time_resolution);
        end
        
        % Selection mask (power percentile + coherence + min speed)
        Pg = sqrt(mP1 .* mP2);                 % symmetric power
        Pthr_local = prctile(Pg(:), meanPowerPercentile);
        candidate_mask = (abs(v_map) >= vmin_mps) & (Pg >= Pthr_local);

        % v1-style blob: threshold on |v|, pick blob with max total power
        CC = bwconncomp(candidate_mask);
        msel = false(NS,NA);
        if CC.NumObjects > 0
            best = 0; ibest = 0;
            for k = 1:CC.NumObjects
                idx = CC.PixelIdxList{k};
                % weight by mean power at t (mP1), like v1’s “strongest blob”
                sP = nansum(mP1(idx));
                if sP > best, best = sP; ibest = k; end
            end
            if ibest > 0
                msel(CC.PixelIdxList{ibest}) = true;
            end
        end
        if ~any(msel(:)), continue; end

        % v1 metric: arithmetic mean of speeds inside the blob (+ within-blob std)
        vvals = v_map(msel);
        zone_frame_speeds{z}(f_idx) = mean(vvals, 'omitnan');
        zone_frame_err{z}(f_idx)    = std(vvals,  'omitnan');

        % per-frame wave-rose for this ROI (optional debug images)
        % if dbgRoiRose_enable && ismember(z, dbgRoiRose_zlist) 
        %     plotRunningRoiRose(z, f_idx, ...
        %         Pg, mP2, mC, phi_int, msel, ...
        %         Scales, Angles, time_resolution, pixel_size_km, ...
        %         roiMask(:,:,z), first_frame_pre, dbgRoiRose_outDir);
        % end
    end
end

% summarize per zone like v1 tracking summary: mean speed and random error
speedMean_mps = nan(Nz,1);
speedStd_mps  = nan(Nz,1);    % temporal random error across frames

for z = 1:Nz
    vseries = zone_frame_speeds{z};
    speedMean_mps(z) = mean(vseries, 'omitnan');
    speedStd_mps(z)  = std(vseries,  'omitnan');  % random error (frame-to-frame)
end

%% === EMPIRICAL SPEED RETRIEVAL WITH FRAME-BY-FRAME ANNOTATION SHIFTS ===

% parameters for 1-D shift (stabilize)
Lwin    = 10;     % max lag searched (bg px)
halfLen = 20;     % half-length of sampled 1-D profile (bg px)

% --- per-zone per-frame outputs (ground & intrinsic)
zone_frame_speed_ground = cell(Nz,1);
zone_frame_speed_intrin = cell(Nz,1);
zone_frame_err_ground   = cell(Nz,1);
zone_frame_err_intrin   = cell(Nz,1);
for z0 = 1:Nz
    zone_frame_speed_ground{z0} = nan(numFrames-1,1);
    zone_frame_speed_intrin{z0} = nan(numFrames-1,1);
    zone_frame_err_ground{z0}   = nan(numFrames-1,1);
    zone_frame_err_intrin{z0}   = nan(numFrames-1,1);
end

px_m = pixel_size_km * 1000;
maxPairs = numFrames - 1;

[LON_sh, LAT_sh] = meshgrid(linspace(lon_vec(1), lon_vec(end), Nx_sh), ...
                            linspace(lat_vec(1), lat_vec(end), Ny_sh));
U_advection_sh = F_U(double(LON_sh), double(LAT_sh));   % [Ny_sh x Nx_sh]
V_advection_sh = F_V(double(LON_sh), double(LAT_sh));   % [Ny_sh x Nx_sh]

for z = 1:Nz
    if isempty(selectedCells{z}) || ~any(selectedCells{z}(:)), continue; end

    % per-Z debug dir
    outDbg = fullfile(wavePacketDir, sprintf('DEBUG_Z%d_LINES', z));
    if ~exist(outDbg,'dir'), mkdir(outDbg); end

    % dims in CWT space (load once outside loops if you prefer)
    S1tmp = load(fullfile(cwtCacheDir, sprintf('cwt_f%06d.mat', 1)), 'W');
    [NyC, NxC, ~, ~] = size(S1tmp.W);

    % ellipse core for this Z
    Ze_full = ellipseMask(zoneProps(z), NyC, NxC);
    Ze_core = Ze_full; % or eroded version
    %Ze_core = imerode(Ze_full, strel('disk', 8));

    % dominant normal for this Z
    [NS, NA] = size(selectedCells{z});
    [~, ai_all] = ind2sub([NS,NA], find(selectedCells{z}));
    th = Angles(ai_all);
    th_dom = atan2(mean(sin(th)), mean(cos(th)));
    khat_en = [cos(th_dom), -sin(th_dom)];
    khat_en = khat_en ./ max(norm(khat_en), eps);
    tang_en = [-khat_en(2), khat_en(1)];
    n = khat_en(:).';

    % per-Z advection fields on CWT and bg grids
    U_advec_cwt = imresize(U_advection_sh, [NyC NxC], 'bilinear');
    V_advec_cwt = imresize(V_advection_sh, [NyC NxC], 'bilinear');

    for f_idx = 1:maxPairs
        % ---- load CWT for the pair
        S1 = load(fullfile(cwtCacheDir, sprintf('cwt_f%06d.mat', f_idx)),   'W');  W1 = S1.W;
        S2 = load(fullfile(cwtCacheDir, sprintf('cwt_f%06d.mat', f_idx+1)), 'W');  W2 = S2.W;
    
        % ---- reconstruct only the selected blob
        A1 = reconstruct_from_cells(W1, selectedCells{z});
        A2 = reconstruct_from_cells(W2, selectedCells{z});
    
        % ---- crest/trough masks inside core
        [crest1, trough1, ok1] = crest_trough_masks(A1, Ze_core, contour_percentile);
        [crest2, trough2, ok2] = crest_trough_masks(A2, Ze_core, contour_percentile);
        %if ~ok1 || ~ok2, continue; end %#ok<*NOPTS>
    
        % ---- inner centroids
        Cc1 = inner_centroids(crest1, Ze_core);
        Ct1 = inner_centroids(trough1, Ze_core);
        Cc2 = inner_centroids(crest2, Ze_core);
        Ct2 = inner_centroids(trough2, Ze_core);

        % ---- rank-locked pairing along normal
        % expected max intrinsic speed (m/s) -> CWT px per frame along n
        vmax = 60;                                % choose per case
        shift_bg = (vmax * time_resolution) / px_m;          % bg px/frame
        shift_cwt = shift_bg * ((khat_en(1)^2)/sx + (khat_en(2)^2)/sy);
        [Pc1,Pc2,idxC] = pair_by_nearest_on_normal(Cc1, Cc2, khat_en, 1.2*shift_cwt);
        [Pt1,Pt2,idxT] = pair_by_nearest_on_normal(Ct1, Ct2, khat_en, 1.2*shift_cwt);
        %if isempty(Pc1) && isempty(Pt1), continue; end
    
        % ---- après avoir obtenu Pc1,Pc2,Pt1,Pt2 (première paire dispo du z)
        Dmed = [median(Pc2-Pc1,1,'omitnan'); median(Pt2-Pt1,1,'omitnan')];
        Dmed = median(Dmed,1,'omitnan');                % vecteur médian en CWT px
        if all(isfinite(Dmed))
            if dot(Dmed, khat_en) < 0                  % khat_en pointe “à contre-sens”
                khat_en = -khat_en;                    % flip une fois pour toutes
                tang_en = [-khat_en(2), khat_en(1)];
            end
        end
    
    
        % ---- backgrounds (uint8) for t and t+Δt
        bg1 = make_preprocessed_bg(fullfile(dataDir, fNames{f_idx}),   varName, ...
                instrument, methodName, fTimes(f_idx), IR_threshold, IR_fillPercentile, ...
                VIS_lowerPercentile, VIS_upperPercentile, VIS_fillPercentile, ...
                clipMinHP, clipMaxHP, lowPassFilterWidth_50, Insolation_Correction, ...
                cumulativeIRmask, doWindow, radius_factor, decay_rate);
        bg2 = make_preprocessed_bg(fullfile(dataDir, fNames{f_idx+1}), varName, ...
                instrument, methodName, fTimes(f_idx+1), IR_threshold, IR_fillPercentile, ...
                VIS_lowerPercentile, VIS_upperPercentile, VIS_fillPercentile, ...
                clipMinHP, clipMaxHP, lowPassFilterWidth_50, Insolation_Correction, ...
                cumulativeIRmask, doWindow, radius_factor, decay_rate);
    
        if shrinkfactor ~= 1; bg1 = imresize(bg1, invshrinkfactor); bg2 = imresize(bg2, invshrinkfactor); end
    
        [NyB, NxB] = size(bg1);
        sx = NxB / NxC;  sy = NyB / NyC;
    
        U_advec_cwt = imresize(U_advection_sh, [NyC NxC], 'bilinear');
        V_advec_cwt = imresize(V_advection_sh, [NyC NxC], 'bilinear');
    
        U_advec_bg = imresize(U_advec_cwt, [NyB NxB], 'bilinear');
        V_advec_bg = imresize(V_advec_cwt, [NyB NxB], 'bilinear');
    
        % ---- centroid-based signed shifts along the common normal (in CWT px)
        % projections centrées sur la NORMALE (toujours dans le bon sens)
        dC = (Pc2 - Pc1) * khat_en(:);   % CWT px
        dT = (Pt2 - Pt1) * khat_en(:);
        
        % robust combine identique
        d_all = [dC; dT];
        %if all(isnan(d_all)), continue; end
        d_med = median(d_all,'omitnan');
        mad   = 1.4826*median(abs(d_all - d_med),'omitnan');
        thr   = max(2, 2.5*mad);
        isInC = ~isnan(dC) & abs(dC - d_med) <= thr;
        isInT = ~isnan(dT) & abs(dT - d_med) <= thr;
        
        % vecteur déplacement médian (CWT px) et vitesses
        dx_px = d_med * khat_en(1);
        dy_px = d_med * khat_en(2);                % rappel: y descend en image
        u_east  =  (dx_px * px_m) / time_resolution;
        v_north = -(dy_px * px_m) / time_resolution;
        v_abs   = hypot(u_east, v_north);       
    
        v_proj_ground = (d_med * px_m) / time_resolution;    % m/s along n
    
        % Per-sample projected shifts (crest+trough), keep inliers only
        d_samples = [dC(isInC); dT(isInT)];                  % CWT px along n
        v_samples_ground = (d_samples * px_m) / time_resolution;  % m/s along n
        
        % Ground speed summary for this frame
        Ubar = mean(U_advec_cwt(Ze_core),'omitnan');
        Vbar = -mean(V_advec_cwt(Ze_core),'omitnan');
        n = khat_en(:).';                     % unit normal you already computed
        U_n = Ubar*n(1) + Vbar*n(2);          % m/s along n
    
        v_samples_intrin = v_samples_ground - U_n;   % sample-wise intrinsic
        
        zone_frame_speed_ground{z}(f_idx) = mean(v_samples_ground,'omitnan');
        zone_frame_err_ground{z}(f_idx)   = std (v_samples_ground,'omitnan');
        zone_frame_speed_intrin{z}(f_idx) = mean(v_samples_intrin,'omitnan');
        zone_frame_err_intrin{z}(f_idx)   = std (v_samples_intrin,'omitnan');
    
        % --- figure with background t, contours at t and t+Δt, centroids, arrows
        fig = figure('Visible','on','Units','pixels','Position',[100 100 NxB NyB]);
        imshow(bg1,[]); hold on;
    
        % --- Mean advection over the core, arrow at core centroid
        S = regionprops(Ze_core, 'Centroid');                 % on CWT grid
        cC = S.Centroid;                                      % [xCWT yCWT]
        x0 = cC(1)*sx;  y0 = cC(2)*sy;                        % bg coords
        
        dx_mean =  Ubar * time_resolution * sx / px_m;        % bg px per frame
        dy_mean =  Vbar * time_resolution * sy / px_m;        % minus for image y
        quiver(x0, y0, dx_mean,  dy_mean, 0, 'c', 'LineWidth', 1.8, 'MaxHeadSize', 2.0);
        
        % --- Along-normal projection of mean wind
        U_n = Ubar*n(1) + Vbar*n(2);                          % m/s along n (geophysical)
        dx_n =  (U_n * time_resolution / px_m) * n(1) * sx;   % bg px per frame
        dy_n =  (U_n * time_resolution / px_m) * n(2) * sy;   % minus for image y
        quiver(x0, y0, dx_n,  dy_n, 0, 'm', 'LineWidth', 1.8, 'MaxHeadSize', 2.0);
    
        % --- Advection field in px/frame on the background grid
        sQ = max(12, round(min(NxB,NyB)/50));
        [Xq,Yq] = meshgrid(1:sQ:NxB, 1:sQ:NyB);
        
        Uq = U_advec_bg(1:sQ:end, 1:sQ:end);   % m/s east
        Vq = - V_advec_bg(1:sQ:end, 1:sQ:end);   % m/s north (geo)
        
        % m/s -> px/frame  (anisotropic pixels; image y va vers le bas)
        Uq_pixpf =  Uq * (time_resolution * sx / px_m);   % Δx en px/frame
        Vq_pixpf =  Vq * (time_resolution * sy / px_m);   % Δy en px/frame
        
        % trace flèches à l’échelle réelle (pas de mise à l’échelle automatique)
        quiver(Xq, Yq, Uq_pixpf,  Vq_pixpf, 0, ...
               'Color',[0.85 0.85 0.85], 'LineWidth',0.8, 'MaxHeadSize',2);
        
        % --- flèche d’échelle: 10 px/frame
        x0=600; y0=40; Lpx=10;
        quiver(x0, y0, Lpx, 0, 0, 'w', 'LineWidth',1.2, 'MaxHeadSize',2);
        text(x0+Lpx+6, y0, '10 px/frame', 'Color','r','FontSize',9, ...
             'VerticalAlignment','middle');    
            
        % Contours at t (solid)
        [B,~]=bwboundaries(crest1,'noholes');
        for k=1:numel(B), b=B{k}; plot(b(:,2)*sx,b(:,1)*sy,'r-','LineWidth',1.2); end
        [B,~]=bwboundaries(trough1,'noholes');
        for k=1:numel(B), b=B{k}; plot(b(:,2)*sx,b(:,1)*sy,'b-','LineWidth',1.2); end
        
        % Contours at t+Δt (dashed)
        [B,~]=bwboundaries(crest2,'noholes');
        for k=1:numel(B), b=B{k}; plot(b(:,2)*sx,b(:,1)*sy,'r--','LineWidth',1.0); end
        [B,~]=bwboundaries(trough2,'noholes');
        for k=1:numel(B), b=B{k}; plot(b(:,2)*sx,b(:,1)*sy,'b--','LineWidth',1.0); end
        
        % Centroids (t) and (t+Δt)
        if ~isempty(Pc1), plot(Pc1(:,1)*sx,Pc1(:,2)*sy,'ro','MarkerSize',6,'LineWidth',1.2); end
        if ~isempty(Pt1), plot(Pt1(:,1)*sx,Pt1(:,2)*sy,'bo','MarkerSize',6,'LineWidth',1.2); end
        if ~isempty(Pc2), plot(Pc2(:,1)*sx,Pc2(:,2)*sy,'rs','MarkerSize',5,'LineWidth',1.0); end
        if ~isempty(Pt2), plot(Pt2(:,1)*sx,Pt2(:,2)*sy,'bs','MarkerSize',5,'LineWidth',1.0); end
        
        % Connectors between t and t+Δt (thin gray)
        if ~isempty(Pc1) && ~isempty(Pc2)
            nConn = min(size(Pc1,1), size(Pc2,1));
            for i=1:nConn
                plot([Pc1(i,1)*sx Pc2(i,1)*sx],[Pc1(i,2)*sy Pc2(i,2)*sy],'-','Color',[0.7 0.7 0.7],'LineWidth',0.8);
            end
        end
        if ~isempty(Pt1) && ~isempty(Pt2)
            nConn = min(size(Pt1,1), size(Pt2,1));
            for i=1:nConn
                plot([Pt1(i,1)*sx Pt2(i,1)*sx],[Pt1(i,2)*sy Pt2(i,2)*sy],'-','Color',[0.7 0.7 0.7],'LineWidth',0.8);
            end
        end
        
        % Arrows: projected shifts along n from t-centroids
        for i=1:numel(dC)
            if ~isfinite(dC(i)) || any(~isfinite(Pc1(i,:))), continue; end
            x0 = Pc1(i,1)*sx; y0 = Pc1(i,2)*sy;
            dx = dC(i)*n(1)*sx; dy = dC(i)*n(2)*sy;
            quiver(x0,y0,dx,dy,0,'Color',[0 1 0],'LineWidth',1.4,'MaxHeadSize',1.2); % green
        end
        for i = 1:numel(dT)
            if ~isfinite(dT(i)) || any(~isfinite(Pt1(i,:))), continue; end
            x0 = Pt1(i,1)*sx; y0 = Pt1(i,2)*sy;
            dx = dT(i)*n(1)*sx; dy = dT(i)*n(2)*sy;
            quiver(x0,y0,dx,dy,0,'Color',[1 .6 0],'LineWidth',1.4,'MaxHeadSize',1.2); % orange
        end
        
        txt = sprintf(['z=%d  f%03d→%03d\n' ...
                       'v_g=%.2f m/s   v_i=%.2f m/s\n' ...
                       'Ū=(%.2f, %.2f) m/s   U_n=%.2f m/s\n' ...
                       'N inliers=%d'], ...
                       z, f_idx, f_idx+1, ...
                       zone_frame_speed_ground{z}(f_idx), ...
                       zone_frame_speed_intrin{z}(f_idx), ...
                       Ubar, Vbar, U_n, ...
                       nnz(isInC)+nnz(isInT));
    
        text(10,20,txt,'Color','w','FontSize',10,'FontWeight','bold', ...
             'BackgroundColor',[0 0 0 0.55],'Margin',5);    

        exportgraphics(fig, fullfile(outDbg, sprintf('Z1_overview_f%03d.png', f_idx)), 'Resolution', 200);
        close(fig);
    
    end
end

%% === SUMMARY FRAME from fixed cells (video module) + v1 stats ===
summaryPath = fullfile(singleOutDir, sprintf('%s_SUMMARY_%s.png', ...
                upper(instrument), datestr(fTimes(1),'yyyymmddHHMM')));

% First-frame background
path_bg = fullfile(dataDir, fNames{1});
raw_bg  = double(ncread(path_bg, varName));
pre_bg  = preprocessFrame(raw_bg, instrument, methodName, ...
                          path_bg, fTimes(1), IR_threshold, IR_fillPercentile, ...
                          VIS_lowerPercentile, VIS_upperPercentile, VIS_fillPercentile, ...
                          clipMinHP, clipMaxHP, lowPassFilterWidth_50, Insolation_Correction, ...
                          cumulativeIRmask);
if doWindow, pre_bg = applyRectangularWindow(pre_bg, radius_factor, decay_rate); end
bg = uint8(255 * mat2gray(pre_bg));  R=bg; G=bg; B=bg;

% CWT of frame #1 to reconstruct contours from fixed cells
S1 = load(fullfile(cwtCacheDir, sprintf('cwt_f%06d.mat', 1)),'W');  W1 = S1.W;
[NyC, NxC, NS, NA] = size(W1);

crest_all  = false(NyC,NxC);
trough_all = false(NyC,NxC);

for z = 1:Nz
    msel = selectedCells{z};
    if isempty(msel) || ~any(msel(:)), continue; end

    % Reconstruct spatial field from the fixed cells (same as video)
    A = zeros(NyC, NxC, 'single');
    [ss, aa] = ind2sub([NS,NA], find(msel));
    for k=1:numel(ss)
        A = A + real(W1(:,:,ss(k),aa(k)));
    end

    % Threshold within the ROI ellipse
    Ze = ellipseMask(zoneProps(z), NyC, NxC);
    vals = A(Ze);
    if isempty(vals) || all(~isfinite(vals)), continue; end
    t = prctile(abs(vals), contour_percentile);

    crest_all  = crest_all  | (Ze & (A >=  t));
    trough_all = trough_all | (Ze & (A <= -t));
end

% Perimeters on top of the background
se = strel('disk', lineW, 0);
crest_perim  = imdilate(bwperim(crest_all),  se);
trough_perim = imdilate(bwperim(trough_all), se);

% Resize perimeters to bg if sizes differ
if ~isequal(size(crest_perim), size(bg))
    crest_perim  = imresize(crest_perim,  size(bg), 'nearest');
    trough_perim = imresize(trough_perim, size(bg), 'nearest');
end

R(crest_perim)  = 255; G(crest_perim)  = 0;   B(crest_perim)  = 0;   % red
R(trough_perim) = 0;   G(trough_perim) = 0;   B(trough_perim) = 255; % blue

fig = figure('Visible','on','Units','pixels','Position',[100 100 size(bg,2) size(bg,1)]);
imshow(cat(3,R,G,B)); hold on;

% --- Scale ROI coordinates from CWT grid to background for labels ---
[NyB, NxB] = size(bg);    
scale_x = NxB / NxC;
scale_y = NyB / NyC;

% Cyan labels with v1-series stats
for z = 1:Nz
    msel = selectedCells{z};
    if isempty(msel) || ~any(msel(:)), continue; end
    if isnan(speedMean_mps(z)) || isnan(speedStd_mps(z)), continue; end
    % scale centroid and axes
    cx = zoneProps(z).Centroid(1)      * scale_x;
    cy = zoneProps(z).Centroid(2)      * scale_y;
    w  = zoneProps(z).MajorAxisLength  * scale_x;
    h  = zoneProps(z).MinorAxisLength  * scale_y;
    anchor = [cx - w/2, max(cy - h/2 - 150, 10)];
    txt = sprintf('Track #%d\nAvg Speed: %.2f m/s\nStd Dev: %.2f m/s', ...
                  z, speedMean_mps(z), speedStd_mps(z));
    text(anchor(1), anchor(2), txt, ...
        'Color',[0 1 1], 'FontWeight','bold', 'FontSize',10, ...
        'HorizontalAlignment','left','VerticalAlignment','top', ...
        'BackgroundColor',[0 0 0 0.55], 'Margin',5);
end
hold off;

exportgraphics(fig, summaryPath, 'Resolution', 300);
close(fig);
fprintf('✓ Summary saved → %s\n', summaryPath);

% Optional cache cleanup
if deleteCWTcacheEnd
    %try rmdir(cwtCacheDir,'s'); fprintf('Deleted CWT cache at %s\n', cwtCacheDir); end
end



%% HELPER FUNCTIONS
%==========================================================================

function drawRoseBlobOutline(ax, maskSA, Scales, Angles)
% Draw red outline of a [NS x NA] logical mask on a rose plotted by DisplayAggregatedWaveRose
    [NS, NA] = size(maskSA);
    if ~any(maskSA(:)), return; end
    per = bwperim(maskSA);
    [s_idx, a_idx] = find(per);
    if isempty(s_idx), return; end

    d_angle = diff(Angles(1:2));
    angle_edges = [Angles - d_angle/2, Angles(end) + d_angle/2];

    log_Scales = log10(Scales);
    d_log = diff(log_Scales(1:2));
    log_edges = [log_Scales(1) - d_log/2, log_Scales + d_log/2];

    hold(ax,'on');
    for k = 1:numel(s_idx)
        si = s_idx(k); ai = a_idx(k);
        theta_c = [angle_edges(ai), angle_edges(ai+1), angle_edges(ai+1), angle_edges(ai)];
        r_c     = [log_edges(si),   log_edges(si),   log_edges(si+1),   log_edges(si+1)];
        [x_c,y_c] = pol2cart(theta_c, r_c);
        patch(ax, x_c, y_c, 'r', 'FaceAlpha', 0, 'EdgeColor', 'r', 'LineWidth', 2);
    end
    hold(ax,'off');
end

function plotRunningRoiRose(zid, f_idx, ...
    meanP1, meanP2, meanC, mean_phi_int, selectedMaskSA, ...
    Scales, Angles, time_resolution, pixel_size_km, ...
    roiMask_z, bg_img, outDir)

    % 1) Convert intrinsic phase (rad) to speed (m/s) for display
    [NS, NA] = size(mean_phi_int);
    mean_speed = zeros(NS,NA,'single');
    for s=1:NS
        lambda_m = Scales(s) * (pixel_size_km*1000) * (pi/sqrt(2));
        mean_speed(s,:) = (mean_phi_int(s,:)/(2*pi)) .* (lambda_m/time_resolution);
    end

    % 2) Figure & panels
    fig = figure('Name',sprintf('ROI %d – Running Rose @ frame %d',zid,f_idx), ...
                 'Visible','off','Units','normalized','OuterPosition',[0 0 1 1]);
    sgtitle(sprintf('ROI %d – Running Averages up to Frame Pair %d',zid, f_idx), ...
            'FontWeight','bold','FontSize',14);

    % Panel 1: Power rose (use P1; you can also show (P1+P2)/2 if you prefer)
    ax1 = subplot(1,3,1);
    DisplayAggregatedWaveRose('Avg Power', meanP1, Scales, Angles, ax1, []);
    title(ax1,'Mean Power Rose');

    % Panel 2: Speed rose (m/s) with blob outline
    ax2 = subplot(1,3,2);
    DisplayAggregatedWaveRose('Speed Intrinsic ', mean_phi_int, Scales, Angles, ax2, [-30, 30]);
    title(ax2,'Mean Intrinsic Speed (m/s)');
    drawRoseBlobOutline(ax2, selectedMaskSA, Scales, Angles);

    % (Optional) annotate a few cells inside the blob
    [si, ai] = ind2sub([NS,NA], find(selectedMaskSA));
    if ~isempty(si)
        maxAnnot = numel(si); 
        lines = cell(1, maxAnnot+1);
        lines{1} = 'Selected blob cells (S, A -> m/s):';
        for k=1:maxAnnot
            s = si(k); a=ai(k);
            lambda_m = Scales(s) * (pixel_size_km*1000) * (pi/sqrt(2));
            v_cell = (mean_phi_int(s,a)/(2*pi)) * (lambda_m/time_resolution);
            lines{k+1} = sprintf('S:%g  A:%.2f  →  %.2f', Scales(s), Angles(a), v_cell);
        end
        pos = get(ax2,'Position');
        annotation('textbox',[pos(1) pos(2)-0.08 pos(3) 0.08], ...
                   'String',lines, 'FontSize',8,'EdgeColor','none', ...
                   'HorizontalAlignment','center','VerticalAlignment','bottom');
    end

    % Panel 3: Spatial ROI overlay on a background image
    ax3 = subplot(1,3,3);
    imagesc(ax3, bg_img); axis(ax3,'image','off'); colormap(ax3,'gray'); hold(ax3,'on');
    red_overlay = cat(3, ones(size(roiMask_z)), zeros(size(roiMask_z)), zeros(size(roiMask_z)));
    h_ov = imshow(red_overlay, 'Parent', ax3);
    set(h_ov,'AlphaData', roiMask_z * 0.4);
    title(ax3, sprintf('ROI %d (running)', zid));
    hold(ax3,'off');

    % Save
    %outfile = fullfile(outDir, sprintf('ROI%03d_runningRose_f%06d.png', zid, f_idx));
    %exportgraphics(fig, outfile, 'Resolution', 200);
    %close(fig);
end

function W = computeCWTframe(img, Scales, Angles)
        w = cwtft2(img, 'wavelet','cauchy', 'scales',Scales,'angles',Angles);
        W = squeeze(w.cfs);
        NS = numel(Scales);
        for s = 1:NS
            norm_factor = (pi/sqrt(2))/Scales(s);
            W(:,:,s,:) = W(:,:,s,:)*norm_factor;
        end
        W = single(W);
    end

function tf = pointInEllipse(pt, E)
    % E must contain: Centroid, MajorAxisLength, MinorAxisLength, Orientation
    cx = E.Centroid(1); cy = E.Centroid(2);
    a  = E.MajorAxisLength/2; b = E.MinorAxisLength/2;
    ang = -E.Orientation * pi/180;  % MATLAB orientation sign convention
    ca = cos(ang); sa = sin(ang);
    x = pt(1) - cx;  y = pt(2) - cy;
    xr =  x*ca + y*sa;
    yr = -x*sa + y*ca;
    tf = (xr.^2)/(a^2) + (yr.^2)/(b^2) <= 1;
end

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

    % --- separate plotting radii from physical scales ---
    %shrinkfactor = evalin('base','shrinkfactor');           % 1 or 2, etc.
    %Scales_plot  = Scales * max(1, shrinkfactor);           % mimics old Scales_orig
    Scales_plot  =  Scales;

    % Log-scale edges (midpoints between log10(Scales))
    log_Scales = log10(Scales_plot);
    d_log      = diff(log_Scales(1:2)); % Assumes uniform log spacing
    log_edges  = [log_Scales(1) - d_log/2, log_Scales + d_log/2];

    % Create meshgrid of edges
    [Theta_edges, R_edges] = meshgrid(angle_edges, log_edges);

   
    if startsWith(labelStr, 'Speed', 'IgnoreCase', true)
        pixel_size_km = evalin('base','pixel_size_km');
        dt            = evalin('base','time_resolution');   % seconds
        for iRow = 1:size(innerpower,1)
            lambda_m = Scales(iRow) * (pixel_size_km*1000) * (pi/sqrt(2));  % wavelet λ
            innerpower(iRow,:) = (innerpower(iRow,:)/(2*pi)) .* (lambda_m/dt); % m/s
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
    max_scale_log = log10(max(Scales_plot));
    for i = 1:length(Scales_plot)
        level_log  = log10(Scales_plot(i));
        theta_ring = linspace(0, 2*pi, 180);
        [x_ring,y_ring] = pol2cart(theta_ring, level_log);
        plot(axHandle, x_ring, y_ring, 'k--', 'LineWidth', 0.5);
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

    % geometry-only scale set for plotting
    shrinkfactor = evalin('base','shrinkfactor');
    Scales_plot  = Scales * max(1, shrinkfactor);

    % Build edges with Scales_plot
    log_Scales = log10(Scales_plot);
    d_log      = diff(log_Scales(1:2));
    log_edges  = [log_Scales(1) - d_log/2, log_Scales + d_log/2];

    % Create meshgrid of edges
    [Theta_edges, R_edges] = meshgrid(angle_edges, log_edges);

    if startsWith(labelStr, 'Speed', 'IgnoreCase', true)
        pixel_size_km = evalin('base','pixel_size_km');
        dt            = evalin('base','time_resolution');   % seconds
        for iRow = 1:size(innerpower,1)
            lambda_m = Scales(iRow) * (pixel_size_km*1000) * (pi/sqrt(2));  % wavelet λ
            innerpower(iRow,:) = (innerpower(iRow,:)/(2*pi)) .* (lambda_m/dt); % m/s
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
    max_scale_log = log10(max(Scales_plot));
    for i = 1:length(Scales_plot)
        level_log  = log10(Scales_plot(i));
        theta_ring = linspace(0-(d_angle/2), pi+(d_angle/2), 180);
        [x_ring,y_ring] = pol2cart(theta_ring, level_log);
        plot(ax1, x_ring, y_ring, 'k--', 'LineWidth', 0.5);
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
        phi_pred(iS,iA) = ((u*dt)/lambda(iS)) * cos(Angles(iA)-alpha) * (pi*2);
    end
end

% Compute the phase difference robustly using complex representation
diffPhase = angle(exp(1i*(pdif - phi_pred)));
% Weight by the coherence (or amplitude)
w = coh;
err = nanmean( (w(:).*diffPhase(:)).^2 );
end

function [roiSpeedCell_corrected, roiResidCell, U_vectors, V_vectors] = applyAdvectionCorrectionROI( ... ...
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
    U_vectors = NaN(numROI, 1);
    V_vectors = NaN(numROI, 1);
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

        % Convert speed and direction to U and V components
        U_vectors(iR) = u_hat * cos(alpha_hat);
        V_vectors(iR) = u_hat * sin(alpha_hat);

        %fprintf('ROI %d: Estimated Advection Speed = %.2f m/s, Direction = %.1f deg\n', ...
        %    iR, u_hat, rad2deg(alpha_hat));

        % 3) Build the theoretical advection phase for each scale, angle
        NSCALES  = numel(Scales);
        NANGLES  = numel(Angles);
        phi_adv  = zeros(NSCALES, NANGLES);
        for iS = 1:NSCALES
            lambda_m = Scales(iS) * (pixel_size_km*1000) * pi/sqrt(2);  % wavelength in meters
            for iA = 1:NANGLES
                phi_adv(iS,iA) = ((u_hat*dt)/lambda_m) * cos(Angles(iA) - alpha_hat) * 2*pi;
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

function drawEllipse(centroid, major_axis, minor_axis, orientation_deg, varargin)
    t=linspace(0,2*pi,100); a=major_axis/2; b=minor_axis/2;
    xy_unit=[a*cos(t);b*sin(t)]; phi=-orientation_deg*pi/180;
    R=[cos(phi),-sin(phi);sin(phi),cos(phi)];
    xy_rot=R*xy_unit;
    plot(xy_rot(1,:)+centroid(1), xy_rot(2,:)+centroid(2), varargin{:});
end

function [fileNames, fileTimestamps, variableName] = getDateRangeFiles(dataDir, startDate, endDate)
    ncFiles=dir(fullfile(dataDir,'*.nc')); if isempty(ncFiles);fileNames={};fileTimestamps=[];variableName='';return; end
    info=ncinfo(fullfile(ncFiles(1).folder,ncFiles(1).name));varList={info.Variables.Name};
    if any(strcmp(varList,'CMI'));variableName='CMI';elseif any(strcmp(varList,'Rad'));variableName='Rad';else;variableName=varList{1};end
    fileNames={};fileTimestamps=datetime([],'ConvertFrom','datenum');
    for i=1:numel(ncFiles);fn=ncFiles(i).name;parts=split(fn,'_');if numel(parts)<6;continue;end
    try;thisDT=datetime(str2double(parts{2}),str2double(parts{3}),str2double(parts{4}),str2double(parts{5}),str2double(erase(parts{6},'.nc')),0);catch;continue;end
    if thisDT>=startDate && thisDT<=endDate;fileNames{end+1}=fn;fileTimestamps(end+1)=thisDT;end;end
    [fileTimestamps,idxSort]=sort(fileTimestamps);fileNames=fileNames(idxSort);
end

function data_preprocessed = preprocessFrame(data, dataType, methodName, fullPath, fileTime, IR_threshold, IR_fillPercentile, VIS_lowerPercentile, VIS_upperPercentile, VIS_fillPercentile, clipMinHP, clipMaxHP, lpWidth50, Insolation_Correction, staticMask)
    switch upper(dataType)
    case 'IR'
        if ~evalin('base','exist(''useCumulativeMask'',''var'')==1') || isempty(staticMask); data(data<IR_threshold)=NaN; nan_mask=isnan(data);
        elseif evalin('base','exist(''useCumulativeMask'',''var'')==1'); nan_mask=staticMask; end
        data(nan_mask)=NaN; data(nan_mask)=mean(data(:),'omitnan');
        data_preprocessed=processDataMethod(data,methodName,clipMinHP,clipMaxHP,lpWidth50);
        fill_value=prctile(data_preprocessed(:),IR_fillPercentile); data_preprocessed(nan_mask')=fill_value;
    otherwise, error('Unsupported data type');
    end
end

function img_processed = processDataMethod(data, methodName, clipMinHP, clipMaxHP, lpWidth50)
    if strcmpi(methodName,'highpass_50_sqrt'); img_p=applyHighPass(data,lpWidth50,true,clipMinHP,clipMaxHP); img_processed=img_p';
    else; error('Unknown method: %s',methodName); end
end

function img_out = applyHighPass(data, filterWidth, doSqrtEnhance, clipMinHP, clipMaxHP)
    lowPass=imgaussfilt(data,filterWidth); highPass=data-lowPass;
    if doSqrtEnhance; highPass=sqrt(abs(highPass)).*sign(highPass); end
    highPass(highPass<clipMinHP)=clipMinHP; highPass(highPass>clipMaxHP)=clipMaxHP;
    img_out=standardizeData(highPass);
end

function data_win = applyRectangularWindow(data_in, radius_factor, decay_rate)
    median_val=median(data_in(:));[rows,cols]=size(data_in);[X,Y]=meshgrid(1:cols,1:rows);
    R=max(abs(X-cols/2)/(cols/2),abs(Y-rows/2)/(rows/2));
    window=1./(1+exp(decay_rate*(R-radius_factor)));
    data_win=window.*data_in+(1-window).*median_val;
end

function out = standardizeData(data)
    if evalin('base','strcmpi(instrument,''IR'') && True_Color_IR'); data=-data; end
    meanVal=mean(data(:),'omitnan');stdVal=std(data(:),'omitnan');out=(data-meanVal)./stdVal;
end

function [angles,projection] = solarPosition(datetime_val,latitude,longitude, time_zone,rotation,dst)
    if iscellstr(datetime_val)||~isvector(datetime_val);datetime_val=datenum(datetime_val);else;datetime_val=datetime_val(:);end
    date=floor(datetime_val);[year,~,~]=datevec(date);time=datetime_val-date;d_n=mod(date-datenum(year,1,1)+1,365);
    B=2*pi*(d_n-81)/365;ET=9.87*sin(2*B)-7.53*cos(B)-1.5*sin(B);
    solarTime=((time*24-double(dst))*60+4*(longitude-time_zone*15)+ET)/60/24;
    latitude_rad=latitude*pi/180;t_h=(solarTime*24-12)*15;t_h_rad=t_h*pi/180;
    delta=-23.45*cos(2*pi*(d_n+10)/365);delta_rad=delta*pi/180;
    theta_rad=acos(sin(latitude_rad)*sin(delta_rad)+cos(latitude_rad)*cos(delta_rad).*cos(t_h_rad));
    theta=theta_rad*180/pi;elevation=90-theta;day=elevation>0;
    cos_phi=(cos(theta_rad)*sin(latitude_rad)-sin(delta_rad))./(sin(theta_rad)*cos(latitude_rad));
    phi_south=acos(min(1,max(-1,cos_phi)));phi_rad=NaN(size(phi_south));phi_rad(day)=pi+sign(t_h(day)).*phi_south(day);
    phi_x=atan2(sin(phi_rad-rotation*pi/180).*sin(theta_rad),cos(theta_rad))*180/pi;
    phi_y=atan2(cos(phi_rad-rotation*pi/180).*sin(theta_rad),cos(theta_rad))*180/pi;
    angles=[theta, phi_rad*180/pi];projection=[phi_x,phi_y];
end

% === IoU between two ellipses (zone vs packet), rasterized ===
function iou = ellipseIoU(Ez, Ep, imsz)
% Ez/Ep: struct with Centroid, MajorAxisLength, MinorAxisLength, Orientation
% imsz: [Ny Nx]
    Ny = imsz(1); Nx = imsz(2);
    mz = ellipseMask(Ez, Ny, Nx);
    mp = ellipseMask(Ep, Ny, Nx);
    inter = nnz(mz & mp);
    uni   = nnz(mz | mp);
    if uni == 0, iou = 0; else, iou = inter / uni; end
end

% === Binary mask from ellipse geometry ===
function m = ellipseMask(E, Ny, Nx)
    [X,Y] = meshgrid(1:Nx, 1:Ny);
    cx = E.Centroid(1); cy = E.Centroid(2);
    a  = E.MajorAxisLength/2; b = E.MinorAxisLength/2;
    ang= -E.Orientation*pi/180;         % MATLAB’s orientation sign
    ca = cos(ang); sa = sin(ang);
    Xc = (X - cx)*ca + (Y - cy)*sa;
    Yc = -(X - cx)*sa + (Y - cy)*ca;
    m  = (Xc.^2)/(a^2) + (Yc.^2)/(b^2) <= 1;
end


function A = reconstruct_from_cells(W, selectedMaskSA)
% W: [Ny x Nx x NS x NA], selectedMaskSA: [NS x NA] logical
    [Ny,Nx,~,~] = size(W);
    A = zeros(Ny,Nx,'single');
    [ss, aa] = ind2sub(size(selectedMaskSA), find(selectedMaskSA));
    for k = 1:numel(ss)
        A = A + real(W(:,:,ss(k),aa(k)));
    end
end

function [crestMask, troughMask, ok] = crest_trough_masks(A, coreMask, percentileVal)
% seuil par percentile de |A| dans coreMask
    v = A(coreMask);
    ok = ~(isempty(v) || all(~isfinite(v)));
    if ~ok
        crestMask = false(size(A));
        troughMask = false(size(A));
        return;
    end
    t = prctile(abs(v), percentileVal);
    crestMask  = coreMask & (A >=  t);
    troughMask = coreMask & (A <= -t);
end

function C = inner_centroids(binaryMask, coreMask)
% centroïdes des composantes connexes non accolées au bord de coreMask
    C = [];
    CC = bwconncomp(binaryMask);
    if CC.NumObjects == 0, return; end
    coreBorder = bwperim(coreMask);
    for q = 1:CC.NumObjects
        idx = CC.PixelIdxList{q};
        %if any(coreBorder(idx)), continue; end
        [yy, xx] = ind2sub(size(binaryMask), idx);
        C = [C; [mean(xx), mean(yy)]]; %#ok<AGROW>
    end
end

function bg = make_preprocessed_bg(ncPath, varName, instrument, methodName, thisTime, ...
    IR_threshold, IR_fillPercentile, VIS_lowerPercentile, VIS_upperPercentile, VIS_fillPercentile, ...
    clipMinHP, clipMaxHP, lowPassFilterWidth_50, Insolation_Correction, cumulativeIRmask, ...
    doWindow, radius_factor, decay_rate)
% produit un fond uint8 [0..255] depuis la frame brute
    raw = double(ncread(ncPath, varName));
    pre = preprocessFrame(raw, instrument, methodName, ...
            ncPath, thisTime, IR_threshold, IR_fillPercentile, ...
            VIS_lowerPercentile, VIS_upperPercentile, VIS_fillPercentile, ...
            clipMinHP, clipMaxHP, lowPassFilterWidth_50, Insolation_Correction, ...
            cumulativeIRmask);
    if doWindow
        pre = applyRectangularWindow(pre, radius_factor, decay_rate);
    end
    bg = uint8(255 * mat2gray(pre));
end

function [P1,P2,idx] = pair_by_nearest_on_normal(C1, C2, k_en, maxGapPx)
% Pair centroids by nearest distance ALONG the common normal.
% - C1: [N1x2], C2: [N2x2] in CWT pixels (x,y)
% - k_en: 1x2 unit normal (east,north) in CWT pixel coordinates
% - maxGapPx (optional): gate in CWT pixels along normal. [] or omit = no gate.
% Returns matched points P1,P2 (same length) and rank index (1..numPairs)
%
% Greedy assignment on 1D projections; robust to births/deaths and reordering.

    P1=[]; P2=[]; idx=[];
    if isempty(C1) || isempty(C2), return; end

    % projection scalars along the normal (CWT px)
    s1 = C1 * k_en(:);
    s2 = C2 * k_en(:);

    % distances along the normal
    D = abs(s1 - s2');                % [N1 x N2]

    % Optional gate
    if nargin < 4 || isempty(maxGapPx), maxGapPx = inf; end
    D(D > maxGapPx) = inf;

    % Greedy nearest matching without crossing (stable enough for 1D)
    used1 = false(size(s1));
    used2 = false(size(s2));
    pairs = [];

    while true
        [dmin, linearIdx] = min(D(:));
        if ~isfinite(dmin), break; end
        [i,j] = ind2sub(size(D), linearIdx);
        pairs(end+1,:) = [i j dmin]; %#ok<AGROW>
        used1(i) = true; used2(j) = true;
        D(i,:) = inf; D(:,j) = inf;   % remove row/col
    end

    if isempty(pairs), return; end

    % Sort pairs by projection rank on frame-1 to keep consistent labeling
    [~,ord] = sort(s1(pairs(:,1)));
    pairs = pairs(ord,:);

    P1  = C1(pairs(:,1),:);
    P2  = C2(pairs(:,2),:);
    idx = (1:size(P1,1)).';
end




