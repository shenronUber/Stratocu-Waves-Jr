% measure_synthetic_wave_error.m (with Advection Removal and improved Peak Detection)
function [perceived_speed, ground_truth_speed, error_ms] = measure_synthetic_wave_error(base_frame, wavelength, direction, cphase, drift_speed_m_s, driftAngleDeg)
%
% This function processes a single synthetic wave case to quantify measurement error.
% It now INCLUDES an advection removal step to isolate the intrinsic wave speed.
%

%% 1) FIXED PARAMETERS FOR THE EXPERIMENT
%----------- SIMULATION SETTINGS ------------------------------------
instrument      = 'IR';
numFrames       = 10;
time_resolution = 1800;

%----------- SYNTHETIC WAVE & PACKET PARAMETERS ---------------------
zamplitude = 100;
PBLdepth   = 1000;
dB_dzPBL   = 0.1;
packet_center_x = 0;
packet_center_y = 0;
packet_width_x  = 4*400e3;
packet_width_y  = 4*300e3;

%----------- SPATIAL SCALING & RESIZING -----------------------------
degrees_per_pixel = 0.04;
km_per_degree     = 111.32;
shrinkfactor      = 2;
invshrinkfactor   = 1 / shrinkfactor;
original_px_km    = degrees_per_pixel * km_per_degree;
pixel_size_km     = original_px_km * shrinkfactor;

%----------- WAVELET PARAMETERS -------------------------------------
Angles = 0 : pi/(7*2) : pi;
Scales = [2, 4, 8, 16, 32, 64];
Scales_orig = Scales;
NSCALES = numel(Scales);
NANGLES = numel(Angles);
if shrinkfactor ~= 1
    Scales = Scales / shrinkfactor;
end

%----------- ROI, PREPROCESSING, & SPEED CORRECTION -----------------
square_size_deg = 10;
methodName      = 'highpass_50_sqrt';
clipMinHP       = -3;
clipMaxHP       = 3;
IR_threshold_setting = 277;
IR_fillPercentile_setting = 50;
useCumulativeMask_setting = false; % We don't use a cumulative mask in this synthetic test
staticMask_data = [];            % Therefore, the static mask is empty
lpWidth50_setting = 50;
beta = 0.8; decayFactor = 0.55; decaySharpness = 1.2; upperCutoff = 16;
nyquistScales= Scales_orig(Scales_orig<=8);
matrixMode = true;
speed_std_factor = 2;
calibrationMatrix = [2, 0.50; 4, 0.50; 8, 0.50; 16, 0.40; 32, 0.335; 64, 0.17; 128, 0.083];

% Parameters for Advection Estimation ---
scalesForAdvection = [4, 8, 16, 32] / shrinkfactor; % Use shrunken scales
cohThreshold = 0.2;

%% 2) GROUND TRUTH CALCULATION (UPDATED)
% Since we are now removing advection, the ground truth we want to recover
% is the wave's own intrinsic phase speed.
ground_truth_speed = cphase;
wave_theta_rad = deg2rad(90 - direction); % Still needed for wave generation and peak selection

%% 3) MAIN PROCESSING LOOP
[rowsF_orig, colsF_orig] = size(base_frame);
rowsF_shrunk = round(rowsF_orig * invshrinkfactor);
colsF_shrunk = round(colsF_orig * invshrinkfactor);

% --- NEW: Add accumulators for coherence calculation ---
crossSpec_sum = zeros(colsF_shrunk, rowsF_shrunk, NSCALES, NANGLES, 'like', 1i);
B1_auto_sum   = zeros(colsF_shrunk, rowsF_shrunk, NSCALES, NANGLES);
B2_auto_sum   = zeros(colsF_shrunk, rowsF_shrunk, NSCALES, NANGLES);
phaseExp_sum  = zeros(colsF_shrunk, rowsF_shrunk, NSCALES, NANGLES, 'like', 1i);
prevWaveletSpec = [];

% Precompute grids for synthetic wave
[X, Y] = meshgrid(1:colsF_orig, 1:rowsF_orig);
DX = 1000 * original_px_km;
Xm = (X - mean(X(:))) * DX;
Ym = (Y - mean(Y(:))) * DX * -1;

for f_idx = 1:numFrames
    data = base_frame;
    k = (2*pi / wavelength) * cos(wave_theta_rad); l = (2*pi / wavelength) * sin(wave_theta_rad); omega = cphase * (2*pi / wavelength); t = (f_idx - 1) * time_resolution; phase = k * Xm + l * Ym - omega * t; dz = zamplitude * sin(phase); Ampwindow = exp(-(((Xm - packet_center_x) / packet_width_x).^2 + ((Ym - packet_center_y) / packet_width_y).^2)); dz = dz .* Ampwindow; dxy = (zamplitude / PBLdepth) * wavelength .* sin(phase - pi/2) ./ DX; dx = dxy .* cos(wave_theta_rad) .* Ampwindow; dy = dxy .* sin(wave_theta_rad) .* Ampwindow; XI = X - dx; YI = Y - dy; warped_img = interp2(X, Y, base_frame, XI, YI, 'linear', 0); data = warped_img .* (1 + dz / PBLdepth * dB_dzPBL);
    if drift_speed_m_s > 0 && f_idx > 1, driftDistance_m = drift_speed_m_s * time_resolution * (f_idx - 1); pxShift = round(driftDistance_m / 1000 / original_px_km); shift_dx = round(pxShift * cosd(driftAngleDeg)); shift_dy = round(-pxShift * sind(driftAngleDeg)); data = circshift(data, [shift_dy, shift_dx]); end
    
    data_pre = preprocessFrame(data, instrument, methodName, IR_threshold_setting, IR_fillPercentile_setting, clipMinHP, clipMaxHP, lpWidth50_setting, useCumulativeMask_setting, staticMask_data);
    data_pre = imresize(data_pre, invshrinkfactor);
    waveStruct = cwtft2(data_pre, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);
    spec_full = squeeze(waveStruct.cfs);
    for iS = 1:NSCALES, spec_full(:,:,iS,:) = spec_full(:,:,iS,:) * ((pi/sqrt(2)) / Scales(iS)); end
    
    % --- Accumulate stats (UPDATED to include coherence terms) ---
    if ~isempty(prevWaveletSpec)
        crossSpec_product = prevWaveletSpec .* conj(spec_full);
        crossSpec_sum = crossSpec_sum + crossSpec_product;
        B1_auto_sum = B1_auto_sum + abs(prevWaveletSpec).^2;
        B2_auto_sum = B2_auto_sum + abs(spec_full).^2;
        phase_mat = angle(crossSpec_product);
        phaseExp_sum = phaseExp_sum + exp(1i * phase_mat);
    end
    prevWaveletSpec = spec_full;
end

%% 4) ROI ANALYSIS & PEAK SPEED EXTRACTION
[rowsF_final, colsF_final] = size(crossSpec_sum, [1 2]);
effective_degrees_per_pixel = degrees_per_pixel * shrinkfactor;
square_size_px = round(square_size_deg / effective_degrees_per_pixel);
center_x = round(rowsF_final / 2); center_y = round(colsF_final / 2);
x_start = max(1, center_x - floor(square_size_px / 2)); y_start = max(1, center_y - floor(square_size_px / 2));
xR = x_start : min(rowsF_final, x_start + square_size_px - 1);
yR = y_start : min(colsF_final, y_start + square_size_px - 1);
numPairs = numFrames - 1;

% Calculate speed-rose for the central ROI
localPhaseExp = phaseExp_sum(xR, yR, :, :);
sumPhaseExp = squeeze(sum(sum(localPhaseExp, 1, 'omitnan'), 2, 'omitnan'));
numPixROI = length(yR) * length(xR);
avgPhase = angle(sumPhaseExp / (numPairs * numPixROI));

% --- NEW: Calculate coherence-rose for the central ROI ---
localCross = crossSpec_sum(xR, yR, :, :);
localB1    = B1_auto_sum(xR, yR, :, :);
localB2    = B2_auto_sum(xR, yR, :, :);
sumCross = squeeze(sum(sum(localCross, 1, 'omitnan'), 2, 'omitnan'));
sumB1    = squeeze(sum(sum(localB1, 1, 'omitnan'), 2, 'omitnan'));
sumB2    = squeeze(sum(sum(localB2, 1, 'omitnan'), 2, 'omitnan'));
gamma_sq = (abs(sumCross).^2) ./ (sumB1 .* sumB2);
roiCohCell = {gamma_sq}; % Put it in a cell for the helper function

% --- NEW: Apply Advection Correction ---
[unadvectedSpeedCell, ~] = applyAdvectionCorrectionROI({avgPhase}, roiCohCell, ...
    Scales, Angles, scalesForAdvection, cohThreshold, pixel_size_km, time_resolution);
unadvected_avgPhase = unadvectedSpeedCell{1};

% --- Correct the UNADFECTED speed-rose and convert to m/s ---
avgPhase_corrected = limitSpeedByScale({unadvected_avgPhase}, Scales_orig, nyquistScales, upperCutoff, beta, decayFactor, decaySharpness, matrixMode, calibrationMatrix);
speedMat_phase = avgPhase_corrected{1};
speedMat_ms = zeros(size(speedMat_phase));
for iRow = 1:NSCALES
    speedMat_ms(iRow, :) = (pixel_size_km * 1000 * (speedMat_phase(iRow, :) / (2*pi)) .* Scales(iRow) * pi/sqrt(2)) / time_resolution;
end

% --- Find the most prominent speed peak (IMPROVED BLOB METHOD) ---
mu = mean(abs(speedMat_ms(:)), 'omitnan'); sigma = std(abs(speedMat_ms(:)), 'omitnan');
threshold = mu + speed_std_factor * sigma;
binary_mask = abs(speedMat_ms) >= threshold;
blobs = bwconncomp(binary_mask);

if blobs.NumObjects == 0
    perceived_speed = NaN;
else
    blob_energies = zeros(blobs.NumObjects, 1);
    for i = 1:blobs.NumObjects, blob_energies(i) = sum(abs(speedMat_ms(blobs.PixelIdxList{i}))); end
    
    [~, best_blob_idx] = max(blob_energies);
    
    % --- NEW: Tie-breaker logic for multiple blobs with similar energy ---
    if blobs.NumObjects > 1
        % Find blobs with energy > 90% of the max energy
        candidate_blob_indices = find(blob_energies > 0.9 * max(blob_energies));
        if numel(candidate_blob_indices) > 1
            % If we have multiple strong candidates, find the one whose
            % average angle is closest to the known wave direction.
            angle_diffs = zeros(size(candidate_blob_indices));
            for i = 1:numel(candidate_blob_indices)
                blob_idx = candidate_blob_indices(i);
                [~, blob_cols] = ind2sub(size(speedMat_ms), blobs.PixelIdxList{blob_idx});
                mean_blob_angle = mean(Angles(blob_cols));
                angle_diffs(i) = abs(mean_blob_angle - wave_theta_rad);
            end
            [~, closest_idx] = min(angle_diffs);
            best_blob_idx = candidate_blob_indices(closest_idx); % Overwrite the winner
        end
    end
    
    winning_blob_indices = blobs.PixelIdxList{best_blob_idx};
    blob_speeds = speedMat_ms(winning_blob_indices);
    [~, peak_in_blob_idx] = max(abs(blob_speeds));
    perceived_speed = blob_speeds(peak_in_blob_idx);
end

% --- Final Error Calculation ---
error_ms = perceived_speed - ground_truth_speed;

end

% =========================================================================
%                         HELPER FUNCTIONS
% =========================================================================

function data_preprocessed = preprocessFrame(data, dataType, methodName, ...
    IR_threshold, IR_fillPercentile, clipMinHP, clipMaxHP, lpWidth50, ...
    useCumulativeMask, staticMask)
% PREPROCESSFRAME (Simplified for synthetic experiment)
%  Applies data-type-specific thresholds and then calls the chosen
%  method-based process.

    switch upper(dataType)
        case 'IR'
            % (1) IR-specific threshold or masking
            if ~useCumulativeMask || isempty(staticMask)
                nan_mask = (data < IR_threshold);
            else
                nan_mask = staticMask;
            end
            data(nan_mask) = NaN;
            
            % Fill masked region with a fallback value
            if any(~nan_mask(:))
                data(nan_mask) = mean(data(~nan_mask), 'omitnan');
            else
                data(:) = 0; % All values are NaN, fill with 0
            end

            % (2) Process by methodName
            data_preprocessed = processDataMethod(data, methodName, ...
                clipMinHP, clipMaxHP, [], lpWidth50, []);

            % (3) Fill masked region in the processed data
            if any(~isnan(data_preprocessed(:)))
                fill_value = prctile(data_preprocessed(:), IR_fillPercentile);
                % The mask needs to be transposed to match the output of processDataMethod
                data_preprocessed(nan_mask') = fill_value;
            end

        otherwise
            error('Unrecognized dataType: %s', dataType);
    end
end
%--------------------------------------------------------------------------

function img_processed = processDataMethod(data, methodName, clipMinHP, clipMaxHP, lpWidth20, lpWidth50, lpWidth100)
    % This function is simplified since we only use 'highpass_50_sqrt'
    switch lower(methodName)
        case 'highpass_50_sqrt'
            img_p = applyHighPass(data, lpWidth50, true, clipMinHP, clipMaxHP);
            img_processed = img_p';
        otherwise
            error('Unknown method: %s', methodName);
    end
end
%--------------------------------------------------------------------------

function img_out = applyHighPass(data, filterWidth, doSqrtEnhance, clipMinHP, clipMaxHP)
    lowPass = imgaussfilt(data, filterWidth, 'FilterDomain', 'spatial');
    highPass = data - lowPass;
    if doSqrtEnhance
        highPass = sqrt(abs(highPass)) .* sign(highPass);
    end
    highPass(highPass < clipMinHP) = clipMinHP;
    highPass(highPass > clipMaxHP) = clipMaxHP;
    img_out = standardizeData(highPass);
end
%--------------------------------------------------------------------------

function out = standardizeData(data)
    meanVal = mean(data(:), 'omitnan');
    stdVal  = std(data(:),  'omitnan');
    if stdVal > 0
        out = (data - meanVal) ./ stdVal;
    else
        out = data - meanVal; % Avoid division by zero if data is flat
    end
end
%--------------------------------------------------------------------------

function roiSpeedCell_final = limitSpeedByScale(roiSpeedCell_in, Scales, nyquistScales, upperCutoff, alpha, decayFactor, decaySharpness, matrixMode, calibrationMatrix)

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
    if any(trustedScalesMask)
        trustedProfile = mean(currentWave(trustedScalesMask, :), 1);
    else
        trustedProfile = zeros(1, nAngles); % Fallback if no trusted scales
    end


    for s_idx = find(nyquistMask')' % Ensure s_idx is a row vector for looping
         s = Scales(s_idx);
         % Detect near-zero phase points
        if ~isempty(find(abs(correctedWave(s_idx,:)) < zeroPhaseThreshold))
                candidates = find(abs(correctedWave(s_idx,:)) < zeroPhaseThreshold);
                if numel(candidates) > 1
                    % Generation of variants
                    variants = cell(numel(candidates), 1);
                    scores = zeros(numel(candidates), 1);

                    for c = 1:numel(candidates)
                        variant = correctedWave(s_idx,:);
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
                    correctedWave(s_idx,:) = variants{bestIdx};

                else

                    [~, anchorIdx] = min(abs(correctedWave(s_idx,:)));  % Finds global minimum index
        
                    % Correct leftward from anchor
                    for a = (anchorIdx-1):-1:1
                        delta = correctedWave(s_idx,a) - correctedWave(s_idx,a+1);
                        if delta > pi
                            correctedWave(s_idx,a) = correctedWave(s_idx,a) - 2*pi;
                        elseif delta < -pi
                            correctedWave(s_idx,a) = correctedWave(s_idx,a) + 2*pi;
                        end
                    end
        
                    % Correct rightward from anchor
                    for a = (anchorIdx+1):nAngles
                        delta = correctedWave(s_idx,a) - correctedWave(s_idx,a-1);
                        if delta > pi
                            correctedWave(s_idx,a) = correctedWave(s_idx,a) - 2*pi;
                        elseif delta < -pi
                            correctedWave(s_idx,a) = correctedWave(s_idx,a) + 2*pi;
                        end
                    end

                end
        else   
            %%% cissors method %%%
            % First pass: left → right
            for a = 2:nAngles
                delta = correctedWave(s_idx,a) - correctedWave(s_idx,a-1);
                if abs(delta) > phaseJumpThreshold
                    if delta > pi
                        correctedWave(s_idx,a) = correctedWave(s_idx,a) - 2*pi;
                    elseif delta < -pi
                        correctedWave(s_idx,a) = correctedWave(s_idx,a) + 2*pi;
                    end
                end
            end
    
            % second pass: right → left (compensation)
            for a = (nAngles-1):-1:1
                delta = correctedWave(s_idx,a) - correctedWave(s_idx,a+1);
                if abs(delta) > phaseJumpThreshold
                    if delta > pi
                        correctedWave(s_idx,a) = correctedWave(s_idx,a) - 2*pi;
                    elseif delta < -pi
                        correctedWave(s_idx,a) = correctedWave(s_idx,a) + 2*pi;
                    end
                end
            end
    
            %%% Result merging %%%
            % Weighted Mean of the two passes 
            correctedWave(s_idx,:) = 0.5*(correctedWave(s_idx,:) + flip(correctedWave(s_idx,:)));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Step 2: Constrained Trend Smoothing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trustedMask = (Scales > max(nyquistScales)) & (Scales <= upperCutoff);
    if any(trustedMask)
        % Calculate baseline trend with scale-dependent decay
        if matrixMode
            baseline = zeros(nScales, nAngles);
            for s_idx = 1:nScales
                thisScale = Scales(s_idx);
                % Interpolate factor from the user-provided table
                factor = interp1(calibrationMatrix(:,1), calibrationMatrix(:,2), thisScale, 'linear', 'extrap');
                baseline(s_idx,:) = mean(correctedWave(trustedMask,:), 1, 'omitnan') * factor;
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
    for s_idx = find(Scales > upperCutoff)'
        correctedWave(s_idx,:) = alpha*baseline(s_idx,:) + (1-alpha)*correctedWave(s_idx,:);
    end

    roiSpeedCell_final{iR} = correctedWave;
end
end

function [roiSpeedCell_corrected, roiResidCell] = applyAdvectionCorrectionROI(roiSpeedCell, roiCohCell, Scales, Angles, scalesForAdvection, cohThreshold, pixel_size_km, dt)
    numROI = numel(roiSpeedCell);
    roiSpeedCell_corrected = cell(size(roiSpeedCell));
    roiResidCell = cell(size(roiSpeedCell));
    for iR = 1:numROI
        pdifrose = roiSpeedCell{iR};
        cohrose = roiCohCell{iR};
        [u_hat, alpha_hat] = estimateAdvectionFromPhase(pdifrose, cohrose, Scales, Angles, scalesForAdvection, cohThreshold, pixel_size_km, dt);
        fprintf('ROI %d: Estimated Advection Speed = %.2f m/s, Direction = %.1f deg\n', iR, u_hat, rad2deg(alpha_hat));
        NSCALES = numel(Scales); NANGLES = numel(Angles);
        phi_adv = zeros(NSCALES, NANGLES);
        for iS = 1:NSCALES
            lambda_m = Scales(iS) * (pixel_size_km*1000) * pi/sqrt(2);
            for iA = 1:NANGLES
                phi_adv(iS,iA) = ((u_hat*dt)/lambda_m) * cos(Angles(iA) - alpha_hat) * pi;
            end
        end
        diffPhase = pdifrose - phi_adv;
        residualPhase = angle( exp(1i * diffPhase) );
        roiSpeedCell_corrected{iR} = residualPhase;
        roiResidCell{iR} = diffPhase;
    end
end

function [u, alpha] = estimateAdvectionFromPhase(pdifrose, cohrose, Scales, Angles, scalesForAdvection, cohThreshold, pixel_size_km, dt)
    mask_scales = ismember(Scales, scalesForAdvection);
    pdif_sel = pdifrose(mask_scales, :);
    coh_sel = sqrt(cohrose(mask_scales, :));
    S_sel = Scales(mask_scales);
    coh_mask = (coh_sel >= cohThreshold);
    pdif_sel(~coh_mask) = NaN;
    lambda_sel = S_sel * (pixel_size_km*1000) * (pi / sqrt(2));
    x0 = [5, pi/2];
    opts = optimset('Display','none');
    x_opt = fminsearch(@(x) costFun(x, pdif_sel, coh_sel, Angles, lambda_sel, dt), x0, opts);
    u = x_opt(1);
    alpha = x_opt(2);
end

function err = costFun(x, pdif, coh, Angles, lambda, dt)
    u = x(1); alpha = x(2);
    [nSc, nAng] = size(pdif);
    phi_pred = zeros(nSc, nAng);
    for iS = 1:nSc
        for iA = 1:nAng
            phi_pred(iS,iA) = ((u*dt)/lambda(iS)) * cos(Angles(iA)-alpha) * pi;
        end
    end
    diffPhase = angle(exp(1i*(pdif - phi_pred)));
    w = coh;
    err = nanmean((w(:).*diffPhase(:)).^2);
end