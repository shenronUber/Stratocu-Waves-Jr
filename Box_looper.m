function SingleFrameWavelet_Code
    %% 1) Basic Setup
    %----------------------------------------------------------------------
    % Date range and folder
    startDate = datetime(2023, 10, 12, 1, 0, 0);
    endDate   = datetime(2023, 10, 12, 9, 0, 0);

    boxRootDir = 'C:\Users\admin\Box\GOES2go_satellite_downloads\4km_SEPAC_IR';  
    waveletResultsDir = 'C:\Users\admin\Box\GWaves_2023_10_11-14_SEPAC\Wavelet_Results'; 

    if ~exist(boxRootDir, 'dir')
        error('Box root directory does not exist: %s', boxRootDir);
    end
    if ~exist(waveletResultsDir, 'dir')
        mkdir(waveletResultsDir);
    end

    % Preprocessing method
    methodName = 'highpass_50_sqrt';  

    % Spatial scaling
    degrees_per_pixel = 0.04;
    km_per_degree     = 111.32;
    original_px_km    = degrees_per_pixel * km_per_degree;  % e.g. 4.4528 km
    shrinkfactor      = 2;  
    invshrinkfactor   = 1 / shrinkfactor;
    pixel_size_km     = original_px_km * shrinkfactor;      % final pixel size after shrink

    % Wavelet parameters
    NANGLES = 12;
    Angles  = 0 : pi/(NANGLES-1) : pi;
    
    Scales_km= [10,50,80,110,140,170,200,230,300,405];
    NSCALES   = numel(Scales_km);
    Scales    = Scales_km / pixel_size_km;   % in pixel units

    % Retrieve all .nc in the date range
    [fileNames, fileTimestamps, variableName] = getDateRangeFiles(boxRootDir, startDate, endDate);

    num_frames = numel(fileNames);
    if num_frames < 1
        error('No files found in the specified range.');
    end
    fprintf('Found %d frames from %s to %s\n', num_frames, ...
        datestr(startDate), datestr(endDate));

    %% 2) Loop Over Frames (Single-Frame Wavelet + Annotated Images)
    %----------------------------------------------------------------------
    for f_idx = 1 : num_frames
        fileName = fileNames{f_idx};
        fileTime = fileTimestamps(f_idx);
        fprintf('[%d/%d] Processing file: %s\n', f_idx, num_frames, fileName);

        % Create subfolder for this frame
        frameDateStr = datestr(fileTime, 'yyyy_mm_dd_HHMMSS');
        outDir = fullfile(waveletResultsDir, sprintf('Frame_%03d_%s', f_idx, frameDateStr));
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end

        % 2a) Read & Preprocess Single Frame
        data = double(ncread(fullfile(boxRootDir, fileName), variableName));
    
        threshold = prctile(data(:), 5);
        data(data < threshold) = NaN;
        nan_mask = isnan(data);
        data(isnan(data)) = mean(data,'all','omitnan');
    
        % Processing chain
        data_pro = processDataMethod(data, methodName);
        data_pro(nan_mask')=prctile(data_pro(:), 50);

        % Apply shrink factor
        if shrinkfactor ~= 1
            data_pro = imresize(data_pro, invshrinkfactor);
        end

        % (Optional) Radial window
        doRadialWindow = true;
        if doRadialWindow
            radius_factor = 1.1;  
            decay_rate    = 0.05; 
            data_pro = applyRadialWindow(data_pro, radius_factor, decay_rate);
        end

        % 2b) Build squares from the final dimension
        % (Equivalent to your "prepare first frame" snippet)
        window_buffer = 10;  % or any user-chosen buffer
        [rowsF, colsF] = size(data_pro);

        x_buffer_range = (window_buffer+1) : (colsF - window_buffer);
        y_buffer_range = (window_buffer+1) : (rowsF - window_buffer);

        adjusted_frame_width  = length(x_buffer_range);
        adjusted_frame_height = length(y_buffer_range);

        % Suppose you want 5° squares. If each pixel is 0.08°, 5° ~ 62.5 px
        % Or you do the version that uses km. For example:
        square_size_deg = 5;  
        effective_deg_per_px = degrees_per_pixel;  % e.g. 0.08
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

        % 2c) Single-Frame Wavelet Transform
        waveStruct = cwtft2(data_pro, 'wavelet','cauchy','scales',Scales,'angles',Angles);
        spec_full  = squeeze(waveStruct.cfs);  % [Ny, Nx, nScales, nAngles]

        % Wavelet Normalisation
        for iS = 1:numel(Scales)
            spec_full(:,:,iS,:) = spec_full(:,:,iS,:) * (2/Scales(iS));
        end

        % Average Wavelet over Scales and Angles
        % Calculate the mean spectrum spatially (averaging over rows and columns)
        spec_full_avg = mean(mean(spec_full, 1), 2);  % Resulting size: [1, 1, nScales, nAngles]
        spec_full_avg= squeeze(spec_full_avg);       % Remove singleton dimensions, size: [nScales, nAngles]

        % Save wavelet .mat
        matName = fullfile(outDir, sprintf('FrameWavelet_%03d.mat', f_idx));
        save(matName, 'spec_full_avg','-v7.3');

        % 2d) Produce (N+1) annotated images: grayscale background + color contours
        produceAnnotatedImages(spec_full, data_pro, squares, Scales, Angles, outDir, f_idx,pixel_size_km,1);
    end

    fprintf('Single-frame wavelet processing completed. Results in %s\n', waveletResultsDir);
end


%% Helper Functions
%==============================================================================
function [fileNames, fileTimestamps, variableName] = getDateRangeFiles(boxDir, startDate, endDate)
    % Gather .nc* files, parse timestamps, filter by date, detect variable

    ncFiles = dir(fullfile(boxDir, '*.nc*'));
    fileNames      = {};
    fileTimestamps = datetime([], 'ConvertFrom','datenum');

    if isempty(ncFiles)
        error('No .nc files in directory: %s', boxDir);
    end

    % Identify variable name from first file
    firstFile = ncFiles(1).name;
    info = ncinfo(fullfile(boxDir, firstFile));
    if any(strcmp({info.Variables.Name}, 'CMI'))
        variableName = 'CMI';
    elseif any(strcmp({info.Variables.Name}, 'Rad'))
        variableName = 'Rad';
    else
        varList = {info.Variables.Name};
        variableName = varList{1};
        warning('No standard var (CMI/Rad). Using %s', variableName);
    end

    % Parse each file
    for i = 1:numel(ncFiles)
        fn = ncFiles(i).name;
        tsStr = extractBetween(fn, '_s', '_e');
        if isempty(tsStr), continue; end
        tsStr = tsStr{1};

        try
            fileTS = datetime(tsStr, 'InputFormat','uuuuDDDHHmmssSSS');
        catch
            continue
        end

        if fileTS >= startDate && fileTS <= endDate
            fileNames{end+1} = fn; 
            fileTimestamps(end+1) = fileTS;
        end
    end

    % Sort
    [fileTimestamps, idxSort] = sort(fileTimestamps);
    fileNames = fileNames(idxSort);
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


function produceAnnotatedImages(spec_full, data_background, squares, ...
                                Scales, Angles, outDir, frameIdx,pixel_size_km,clevfactor)
    % Produce (N+1) images: one per scale, plus one "summed" image
    [Ny, Nx, nScales, nAngles] = size(spec_full);
    power = abs(spec_full) .^2;
    buffer = round( max(Scales) );
    innerpower = squeeze( mean(mean( power(buffer:size(power,1)-buffer, ...
                                   buffer:size(power,2)-buffer, :,:) )));
        

    % (1) For each scale
    for s_idx = 1 : nScales

        % Figure: show original data in grayscale, overlay wave_sumAngles as contours
        fig = figure('visible','off');
        
        % Show background in grayscale
        imagesc(data_background, [0, 1]);  % or auto colormap
        colormap(gray);
        axis image off;
        hold on;

        % Trouver l'angle correspondant à la valeur maximale de power
        [max_value, max_idx_angle] = max(innerpower(s_idx,:));  

        %if max_value < mean(innerpower(s_idx,:))*1.5
        %   continue
        %end

        % Extract the wavelet coefficients at the specified scale and angle
        wavelet_real = real(spec_full(:, :, s_idx, max_idx_angle));
        wavelet_abs = abs(spec_full(:, :, s_idx, max_idx_angle));
    
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

        % Overlay squares
        for sq = 1 : numel(squares)
            rectangle('Position',...
                [squares(sq).x_range(1), squares(sq).y_range(1), ...
                 length(squares(sq).x_range), length(squares(sq).y_range)], ...
                 'EdgeColor','yellow','LineWidth',1);
        end

        % Title with frame number & scale in km
        wavelength_km = Scales(s_idx)*pi/sqrt(2)*pixel_size_km;
        titleStr = sprintf('Frame %d, Wavelength ~ %.1f km', frameIdx, wavelength_km );
        title({titleStr}, 'Color','k','FontWeight','bold','FontSize',20);
        subtitleStr=sprintf('Angle : %f °',Angles(max_idx_angle)*180/pi);
        subtitle({subtitleStr}, 'Color','k','FontSize',10);

        outName = fullfile(outDir, sprintf('Frame%03d_Scale%02d.png', frameIdx, s_idx));
        saveas(fig, outName);
        close(fig);
    end
end
