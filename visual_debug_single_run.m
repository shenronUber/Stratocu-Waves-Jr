% visual_debug_single_run.m (ADAPTED FOR ADVECTION REMOVAL)
%
% This script runs a single synthetic wave test and generates visual outputs
% to debug the FULL pipeline, including advection estimation and removal.

clear; clc; close all;

%% 1. DEFINE THE TEST CASE TO VISUALLY DEBUG
test_params.wavelength      = 150e3;    % meters
test_params.direction       = 235;      % degrees
test_params.cphase          = 15;       % m/s
test_params.drift_speed_m_s = 10;       % m/s
test_params.driftAngleDeg   = 90;       % degrees

% --- Setup the Environment ---
base_frame_path = 'C:\Users\admin\Box\GWaves_2023_10_11-14_SEPAC\IR\Data\IR_2023_10_12_14_15.nc'; % <<< EDIT THIS PATH
varName         = 'CMI';

% --- Create an output directory for the debug images ---
outputDir = 'debug_output_adv_removal';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
fprintf('Debug outputs will be saved in: ./%s/\n', outputDir);

%% 2. LOAD BASE FRAME
base_frame = double(ncread(base_frame_path, varName));
fprintf('Loaded base frame from: %s\n', base_frame_path);

%% 3. RUN THE ANALYSIS (Embedded logic from the LATEST measure_synthetic_wave_error)

% --- Start of embedded logic ---
% Fixed parameters
instrument = 'IR'; numFrames = 10; time_resolution = 1800; zamplitude = 100; PBLdepth = 1000; dB_dzPBL = 0.1; packet_center_x = 0; packet_center_y = 0; packet_width_x = 4*400e3; packet_width_y = 4*300e3; degrees_per_pixel = 0.04; km_per_degree = 111.32; shrinkfactor = 2; invshrinkfactor = 1 / shrinkfactor; original_px_km = degrees_per_pixel * km_per_degree; pixel_size_km = original_px_km * shrinkfactor; Angles = 0 : pi/(7*2) : pi; Scales = [2, 4, 8, 16, 32, 64]; Scales_orig = Scales; NSCALES = numel(Scales); NANGLES = numel(Angles); if shrinkfactor ~= 1, Scales = Scales / shrinkfactor; end; square_size_deg = 10; methodName = 'highpass_50_sqrt'; clipMinHP = -3; clipMaxHP = 3; IR_threshold_setting = 277; IR_fillPercentile_setting = 50; useCumulativeMask_setting = false; staticMask_data = []; lpWidth50_setting = 50; beta = 0.8; decayFactor = 0.55; decaySharpness = 1.2; upperCutoff = 16; nyquistScales= Scales_orig(Scales_orig<=8); matrixMode = true; speed_std_factor = 2; calibrationMatrix = [2,0.5; 4,0.5; 8,0.5; 16,0.4; 32,0.335; 64,0.17; 128,0.083];
scalesForAdvection = [4, 8, 16, 32] / shrinkfactor; cohThreshold = 0.2;

% UPDATED Ground Truth Calculation
ground_truth_speed = test_params.cphase;
wave_theta_rad = deg2rad(90 - test_params.direction);

% Main Processing Loop
[rowsF_orig, colsF_orig] = size(base_frame); rowsF_shrunk = round(rowsF_orig * invshrinkfactor); colsF_shrunk = round(colsF_orig * invshrinkfactor);
crossSpec_sum = zeros(colsF_shrunk, rowsF_shrunk, NSCALES, NANGLES, 'like', 1i); B1_auto_sum = zeros(colsF_shrunk, rowsF_shrunk, NSCALES, NANGLES); B2_auto_sum = zeros(colsF_shrunk, rowsF_shrunk, NSCALES, NANGLES); phaseExp_sum = zeros(colsF_shrunk, rowsF_shrunk, NSCALES, NANGLES, 'like', 1i); prevWaveletSpec = [];
[X, Y] = meshgrid(1:colsF_orig, 1:rowsF_orig); DX = 1000 * original_px_km; Xm = (X - mean(X(:))) * DX; Ym = (Y - mean(Y(:))) * DX * -1;
for f_idx = 1:numFrames
    data = base_frame; k = (2*pi / test_params.wavelength) * cos(wave_theta_rad); l = (2*pi / test_params.wavelength) * sin(wave_theta_rad); omega = test_params.cphase * (2*pi / test_params.wavelength); t = (f_idx - 1) * time_resolution; phase = k * Xm + l * Ym - omega * t; dz = zamplitude * sin(phase); Ampwindow = exp(-(((Xm - packet_center_x) / packet_width_x).^2 + ((Ym - packet_center_y) / packet_width_y).^2)); dz = dz .* Ampwindow; dxy = (zamplitude / PBLdepth) * test_params.wavelength .* sin(phase - pi/2) ./ DX; dx = dxy .* cos(wave_theta_rad) .* Ampwindow; dy = dxy .* sin(wave_theta_rad) .* Ampwindow; XI = X - dx; YI = Y - dy; warped_img = interp2(X, Y, base_frame, XI, YI, 'linear', 0); data = warped_img .* (1 + dz / PBLdepth * dB_dzPBL);
    if test_params.drift_speed_m_s > 0 && f_idx > 1, driftDistance_m = test_params.drift_speed_m_s * time_resolution * (f_idx - 1); pxShift = round(driftDistance_m / 1000 / original_px_km); shift_dx = round(pxShift * cosd(test_params.driftAngleDeg)); shift_dy = round(-pxShift * sind(test_params.driftAngleDeg)); data = circshift(data, [shift_dy, shift_dx]); end
    if f_idx == 1, figure('Visible', 'off'); imagesc(data'); colormap(gray); axis image off; title('Generated Frame 1'); saveas(gcf, fullfile(outputDir, '01_generated_frame_01.png')); end
    if f_idx == numFrames, figure('Visible', 'off'); imagesc(data'); colormap(gray); axis image off; title(sprintf('Generated Frame %d', numFrames)); saveas(gcf, fullfile(outputDir, '02_generated_frame_last.png'));
        data_pre_final_vis = preprocessFrame(data, instrument, methodName, IR_threshold_setting, IR_fillPercentile_setting, clipMinHP, clipMaxHP, lpWidth50_setting, useCumulativeMask_setting, staticMask_data);
        figure('Visible', 'off'); imagesc(data_pre_final_vis'); colormap(gray); axis image off; title('Final Pre-Processed Frame'); saveas(gcf, fullfile(outputDir, '03_preprocessed_frame_last.png'));
    end
    data_pre = preprocessFrame(data, instrument, methodName, IR_threshold_setting, IR_fillPercentile_setting, clipMinHP, clipMaxHP, lpWidth50_setting, useCumulativeMask_setting, staticMask_data);
    data_pre = imresize(data_pre, invshrinkfactor);
    waveStruct = cwtft2(data_pre, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles); spec_full = squeeze(waveStruct.cfs);
    for iS = 1:NSCALES, spec_full(:,:,iS,:) = spec_full(:,:,iS,:) * ((pi/sqrt(2)) / Scales(iS)); end
    if ~isempty(prevWaveletSpec), crossSpec_product = prevWaveletSpec .* conj(spec_full); crossSpec_sum = crossSpec_sum + crossSpec_product; B1_auto_sum = B1_auto_sum + abs(prevWaveletSpec).^2; B2_auto_sum = B2_auto_sum + abs(spec_full).^2; phase_mat = angle(crossSpec_product); phaseExp_sum = phaseExp_sum + exp(1i * phase_mat); end
    prevWaveletSpec = spec_full;
end

% ROI Analysis
[rowsF_final, colsF_final] = size(crossSpec_sum, [1 2]); effective_degrees_per_pixel = degrees_per_pixel * shrinkfactor; square_size_px = round(square_size_deg / effective_degrees_per_pixel); center_x = round(rowsF_final / 2); center_y = round(colsF_final / 2); x_start = max(1, center_x - floor(square_size_px / 2)); y_start = max(1, center_y - floor(square_size_px / 2));
xR = x_start : min(rowsF_final, x_start + square_size_px - 1); yR = y_start : min(colsF_final, y_start + square_size_px - 1);
figure('Visible', 'off'); imagesc(data_pre); colormap(gray); axis image off; hold on; rectangle('Position', [yR(1), xR(1), numel(yR), numel(xR)], 'EdgeColor', 'r', 'LineWidth', 2); title('Central ROI for Analysis'); saveas(gcf, fullfile(outputDir, '04_central_roi.png'));
numPairs = numFrames - 1;

% Total Speed Rose (Before Advection Removal)
localPhaseExp = phaseExp_sum(xR, yR, :, :); sumPhaseExp = squeeze(sum(sum(localPhaseExp, 1, 'omitnan'), 2, 'omitnan')); numPixROI = length(yR) * length(xR);
avgPhase_total = angle(sumPhaseExp / (numPairs * numPixROI));

% --- VISUALIZATION POINT 5: Plot Total Speed-Rose ---
avgPhase_total_corr = limitSpeedByScale({avgPhase_total}, Scales_orig, nyquistScales, upperCutoff, beta, decayFactor, decaySharpness, matrixMode, calibrationMatrix);
speedMat_total_ms = zeros(size(avgPhase_total_corr{1}));
for iRow = 1:NSCALES, speedMat_total_ms(iRow, :) = (pixel_size_km * 1000 * (avgPhase_total_corr{1}(iRow, :) / (2*pi)) .* Scales(iRow) * pi/sqrt(2)) / time_resolution; end
figure('Visible', 'off'); [Theta, R] = meshgrid(Angles, log10(Scales_orig)); [X, Y] = pol2cart(Theta, R); pcolor(X, Y, speedMat_total_ms); shading interp; colormap('jet'); colorbar; axis equal tight off;
title('Total Speed-Rose (Before Advection Removal)');
saveas(gcf, fullfile(outputDir, '05_total_speed_rose.png'));

% Coherence Rose Calculation
localCross = crossSpec_sum(xR, yR, :, :); localB1 = B1_auto_sum(xR, yR, :, :); localB2 = B2_auto_sum(xR, yR, :, :);
sumCross = squeeze(sum(sum(localCross, 1, 'omitnan'), 2, 'omitnan')); sumB1 = squeeze(sum(sum(localB1, 1, 'omitnan'), 2, 'omitnan')); sumB2 = squeeze(sum(sum(localB2, 1, 'omitnan'), 2, 'omitnan'));
gamma_sq = (abs(sumCross).^2) ./ (sumB1 .* sumB2);
roiCohCell = {gamma_sq};

% Advection Correction
[unadvectedSpeedCell, ~] = applyAdvectionCorrectionROI({avgPhase_total}, roiCohCell, Scales, Angles, scalesForAdvection, cohThreshold, pixel_size_km, time_resolution);
unadvected_avgPhase = unadvectedSpeedCell{1};

% Unadvected (Final) Speed-Rose
avgPhase_corrected = limitSpeedByScale({unadvected_avgPhase}, Scales_orig, nyquistScales, upperCutoff, beta, decayFactor, decaySharpness, matrixMode, calibrationMatrix);
speedMat_phase = avgPhase_corrected{1};
speedMat_ms = zeros(size(speedMat_phase));
for iRow = 1:NSCALES, speedMat_ms(iRow, :) = (pixel_size_km * 1000 * (speedMat_phase(iRow, :) / (2*pi)) .* Scales(iRow) * pi/sqrt(2)) / time_resolution; end

% Final Peak Detection
mu = mean(abs(speedMat_ms(:)), 'omitnan'); sigma = std(abs(speedMat_ms(:)), 'omitnan'); threshold = mu + speed_std_factor * sigma; binary_mask = abs(speedMat_ms) >= threshold; blobs = bwconncomp(binary_mask);
if blobs.NumObjects == 0, perceived_speed = NaN; max_idx = NaN;
else
    blob_energies = zeros(blobs.NumObjects, 1); for i = 1:blobs.NumObjects, blob_energies(i) = sum(abs(speedMat_ms(blobs.PixelIdxList{i}))); end
    [~, best_blob_idx] = max(blob_energies);
    if blobs.NumObjects > 1, candidate_blob_indices = find(blob_energies > 0.9 * max(blob_energies));
        if numel(candidate_blob_indices) > 1, angle_diffs = zeros(size(candidate_blob_indices));
            for i = 1:numel(candidate_blob_indices)
                blob_idx = candidate_blob_indices(i); [~, blob_cols] = ind2sub(size(speedMat_ms), blobs.PixelIdxList{blob_idx});
                mean_blob_angle = mean(Angles(blob_cols)); angle_diffs(i) = abs(mean_blob_angle - wave_theta_rad);
            end
            [~, closest_idx] = min(angle_diffs); best_blob_idx = candidate_blob_indices(closest_idx);
        end
    end
    winning_blob_indices = blobs.PixelIdxList{best_blob_idx}; blob_speeds = speedMat_ms(winning_blob_indices); [~, peak_in_blob_idx] = max(abs(blob_speeds));
    perceived_speed = blob_speeds(peak_in_blob_idx); max_idx = winning_blob_indices(peak_in_blob_idx);
end
error_ms = perceived_speed - ground_truth_speed;

% --- VISUALIZATION POINT 6: Plot Final Unadvected Speed-Rose ---
figure('Visible', 'off'); [Theta, R] = meshgrid(Angles, log10(Scales_orig)); [X, Y] = pol2cart(Theta, R); pcolor(X, Y, speedMat_ms); shading interp; colormap('jet'); colorbar; axis equal tight off; hold on;
if ~isnan(perceived_speed), [peak_row, peak_col] = ind2sub(size(speedMat_ms), max_idx); peak_angle = Angles(peak_col); peak_scale = log10(Scales_orig(peak_row)); [peak_x, peak_y] = pol2cart(peak_angle, peak_scale); plot(peak_x, peak_y, 'kx', 'MarkerSize', 15, 'LineWidth', 3); plot(peak_x, peak_y, 'wo', 'MarkerSize', 15, 'LineWidth', 2); end
title(sprintf('Final (Unadvected) Speed-Rose (m/s)\nDetected Peak: %.2f m/s', perceived_speed));
saveas(gcf, fullfile(outputDir, '06_final_speed_rose_with_peak.png'));

close all;
fprintf('--- Debug complete. ---\n');
fprintf('Ground Truth (intrinsic wave speed): %.2f m/s\n', ground_truth_speed);
fprintf('Perceived Speed (after advection removal): %.2f m/s\n', perceived_speed);
fprintf('Error: %.2f m/s\n', error_ms);


% =========================================================================
%                         HELPER FUNCTIONS
% =========================================================================


function data_preprocessed = preprocessFrame(data, dataType, methodName, ...
    IR_threshold, IR_fillPercentile, clipMinHP, clipMaxHP, lpWidth50, ...
    useCumulativeMask, staticMask)
% ... (full function code)
    switch upper(dataType)
        case 'IR'
            if ~useCumulativeMask || isempty(staticMask)
                nan_mask = (data < IR_threshold);
            else
                nan_mask = staticMask;
            end
            data(nan_mask) = NaN;
            if any(~nan_mask(:)), data(nan_mask) = mean(data(~nan_mask), 'omitnan'); else, data(:) = 0; end
            data_preprocessed = processDataMethod(data, methodName, clipMinHP, clipMaxHP, [], lpWidth50, []);
            if any(~isnan(data_preprocessed(:)))
                fill_value = prctile(data_preprocessed(:), IR_fillPercentile);
                data_preprocessed(nan_mask') = fill_value;
            end
        otherwise, error('Unrecognized dataType: %s', dataType);
    end
end
function img_processed = processDataMethod(data, methodName, clipMinHP, clipMaxHP, lpWidth20, lpWidth50, lpWidth100)
    switch lower(methodName)
        case 'highpass_50_sqrt'
            img_p = applyHighPass(data, lpWidth50, true, clipMinHP, clipMaxHP);
            img_processed = img_p';
        otherwise, error('Unknown method: %s', methodName);
    end
end
function img_out = applyHighPass(data, filterWidth, doSqrtEnhance, clipMinHP, clipMaxHP)
    lowPass = imgaussfilt(data, filterWidth, 'FilterDomain', 'spatial');
    highPass = data - lowPass;
    if doSqrtEnhance, highPass = sqrt(abs(highPass)) .* sign(highPass); end
    highPass(highPass < clipMinHP) = clipMinHP;
    highPass(highPass > clipMaxHP) = clipMaxHP;
    img_out = standardizeData(highPass);
end
function out = standardizeData(data)
    meanVal = mean(data(:), 'omitnan'); stdVal = std(data(:), 'omitnan');
    if stdVal > 0, out = (data - meanVal) ./ stdVal; else, out = data - meanVal; end
end
function roiSpeedCell_final = limitSpeedByScale(roiSpeedCell_in, Scales, nyquistScales, upperCutoff, alpha, decayFactor, decaySharpness, matrixMode, calibrationMatrix)
% ... (PASTE THE FULL, CORRECT limitSpeedByScale function here) ...
nROI = numel(roiSpeedCell_in);
roiSpeedCell_final = cell(size(roiSpeedCell_in));
phaseJumpThreshold = pi * 0.6;
zeroPhaseThreshold = 1;
for iR = 1:nROI
    currentWave = roiSpeedCell_in{iR};
    [nScales, nAngles] = size(currentWave);
    nyquistMask = ismember(Scales, nyquistScales);    
    correctedWave = currentWave;
    trustedScalesMask = ~ismember(Scales, nyquistScales);
    if any(trustedScalesMask), trustedProfile = mean(currentWave(trustedScalesMask, :), 1);
    else, trustedProfile = zeros(1, nAngles); end
    for s_idx = find(nyquistMask')'
         if ~isempty(find(abs(correctedWave(s_idx,:)) < zeroPhaseThreshold))
                candidates = find(abs(correctedWave(s_idx,:)) < zeroPhaseThreshold);
                if numel(candidates) > 1
                    variants = cell(numel(candidates), 1); scores = zeros(numel(candidates), 1);
                    for c = 1:numel(candidates)
                        variant = correctedWave(s_idx,:); anchor = candidates(c);
                        for a = anchor-1:-1:1
                            delta = variant(a+1) - variant(a);
                            if abs(delta) > phaseJumpThreshold, variant(a) = variant(a) + round(delta/(2*pi))*2*pi; end
                        end
                        for a = anchor+1:nAngles
                            delta = variant(a) - variant(a-1);
                            if abs(delta) > phaseJumpThreshold, variant(a) = variant(a) - round(delta/(2*pi))*2*pi; end
                        end
                        scores(c) = sqrt(mean((variant - trustedProfile).^2)); variants{c} = variant;
                    end
                    [~, bestIdx] = min(scores); correctedWave(s_idx,:) = variants{bestIdx};
                else
                    [~, anchorIdx] = min(abs(correctedWave(s_idx,:)));
                    for a = (anchorIdx-1):-1:1
                        delta = correctedWave(s_idx,a) - correctedWave(s_idx,a+1);
                        if delta > pi, correctedWave(s_idx,a) = correctedWave(s_idx,a) - 2*pi;
                        elseif delta < -pi, correctedWave(s_idx,a) = correctedWave(s_idx,a) + 2*pi; end
                    end
                    for a = (anchorIdx+1):nAngles
                        delta = correctedWave(s_idx,a) - correctedWave(s_idx,a-1);
                        if delta > pi, correctedWave(s_idx,a) = correctedWave(s_idx,a) - 2*pi;
                        elseif delta < -pi, correctedWave(s_idx,a) = correctedWave(s_idx,a) + 2*pi; end
                    end
                end
        else
            for a = 2:nAngles
                delta = correctedWave(s_idx,a) - correctedWave(s_idx,a-1);
                if abs(delta) > phaseJumpThreshold
                    if delta > pi, correctedWave(s_idx,a) = correctedWave(s_idx,a) - 2*pi;
                    elseif delta < -pi, correctedWave(s_idx,a) = correctedWave(s_idx,a) + 2*pi; end
                end
            end
            for a = (nAngles-1):-1:1
                delta = correctedWave(s_idx,a) - correctedWave(s_idx,a+1);
                if abs(delta) > phaseJumpThreshold
                    if delta > pi, correctedWave(s_idx,a) = correctedWave(s_idx,a) - 2*pi;
                    elseif delta < -pi, correctedWave(s_idx,a) = correctedWave(s_idx,a) + 2*pi; end
                end
            end
            correctedWave(s_idx,:) = 0.5*(correctedWave(s_idx,:) + flip(correctedWave(s_idx,:)));
        end
    end
    trustedMask = (Scales > max(nyquistScales)) & (Scales <= upperCutoff);
    if any(trustedMask)
        if matrixMode
            baseline = zeros(nScales, nAngles);
            for s_idx = 1:nScales
                thisScale = Scales(s_idx);
                factor = interp1(calibrationMatrix(:,1), calibrationMatrix(:,2), thisScale, 'linear', 'extrap');
                baseline(s_idx,:) = mean(correctedWave(trustedMask,:), 1, 'omitnan') * factor;
            end
        else
            baseline = mean(correctedWave(trustedMask,:), 1, 'omitnan') .* exp(-decaySharpness * (Scales'/max(Scales))) * decayFactor;
        end
    else, baseline = zeros(nScales, 1); end
    for s_idx = find(Scales > upperCutoff)', correctedWave(s_idx,:) = alpha*baseline(s_idx,:) + (1-alpha)*correctedWave(s_idx,:); end
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