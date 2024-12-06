%% Parameters and Setup

% Define parameters
startDate = datetime(2023, 10, 11);
endDate = datetime(2023, 10, 14);
downloadDir = 'downloaded_data/IR';
dataType = 'IR';

% Define methods you want to produce videos for
methods = {
    'raw_normalized',    'IR_20231011_20231014_raw_normalized.mj2';
    'truncated',         'IR_20231011_20231014_truncated.mj2';
    'highpass_20',       'IR_20231011_20231014_highpass_20.mj2';
    'highpass_20_sqrt',  'IR_20231011_20231014_highpass_20_sqrt.mj2';
    'highpass_50',       'IR_20231011_20231014_highpass_50.mj2';
    'highpass_50_sqrt',  'IR_20231011_20231014_highpass_50_sqrt.mj2';
    'highpass_100',      'IR_20231011_20231014_highpass_100.mj2';
    'highpass_100_sqrt', 'IR_20231011_20231014_highpass_100_sqrt.mj2';
    'raw_masked',        'IR_20231011_20231014_raw_masked.mj2';
};

%% List and Select Files Within Date Range

% List all .nc files in the download directory
ncFiles = dir(fullfile(downloadDir, '*.nc*'));  % Includes .nc and .nc.4km files

if isempty(ncFiles)
    error('No .nc files found in the download directory.');
end

% Initialize arrays to store filenames and timestamps
fileNames = {};
fileTimestamps = datetime([], 'ConvertFrom', 'datenum');

for i = 1:length(ncFiles)
    fileName = ncFiles(i).name;
    
    % Extract timestamp between '_s' and '_e'
    timestampStr = extractBetween(fileName, '_s', '_e');
    if isempty(timestampStr)
        continue;
    end
    
    timestampStr = timestampStr{1};
    % The timestamp format is 'yyyyDDDHHmmssSSS'
    timestampFormat = 'uuuuDDDHHmmssSSS'; 
    
    try
        fileTimestamp = datetime(timestampStr, 'InputFormat', timestampFormat);
        % Check if the file timestamp is within the specified date range
        if fileTimestamp >= startDate && fileTimestamp <= endDate
            fileNames{end+1} = fileName;
            fileTimestamps(end+1) = fileTimestamp;
        end
    catch
        continue;
    end
end

% Check if any files were found within the date range
if isempty(fileNames)
    error('No files found within the specified date range.');
end

% Sort files based on extracted timestamps
[fileTimestamps, sortIdx] = sort(fileTimestamps);
fileNames = fileNames(sortIdx);

% Determine the variable name (CMI or Rad) from the first file
filePath = fullfile(downloadDir, fileNames{1});
info = ncinfo(filePath);
if any(strcmp({info.Variables.Name}, 'CMI'))
    variableName = 'CMI';
elseif any(strcmp({info.Variables.Name}, 'Rad'))
    variableName = 'Rad';
else
    availableVars = {info.Variables.Name};
    warning('No suitable variable found. Available variables: %s', strjoin(availableVars, ', '));
end

%% Preview Lab: Display the first frame for each method

% Load the first frame's data
firstFilePath = fullfile(downloadDir, fileNames{1});
data = double(ncread(firstFilePath, variableName));

for m = 1:size(methods, 1)
    methodName = methods{m,1};
    % Process single frame for preview
    img_preview = processDataMethod(data, methodName);
    figure('Name', sprintf('Preview - %s', methodName));
    imshow(img_preview)
    title(methodName, 'Interpreter', 'none');
end

%% After reviewing the previews, create videos for each method
% (You can comment out the video creation step until satisfied with the preview)

for m = 1:size(methods, 1)
    methodName = methods{m,1};
    outputVideoFile = methods{m,2};
    fprintf('Creating video for method: %s -> %s\n', methodName, outputVideoFile);
    createVideoForMethod(fileNames, downloadDir, dataType, variableName, outputVideoFile, methodName);
end

fprintf('All videos created.\n');

%% Function Definitions

function createVideoForMethod(fileNames, downloadDir, dataType, variableName, outputVideoFile, methodName)
    v = VideoWriter(outputVideoFile, 'Motion JPEG 2000');
    v.FrameRate = 2;  % Adjust frame rate as needed
    v.LosslessCompression = true; % Enable lossless compression
    open(v);

    for i = 1:length(fileNames)
        fileName = fileNames{i};
        filePath = fullfile(downloadDir, fileName);
        fprintf('Processing %s with method %s...\n', fileName, methodName);

        try
            data = double(ncread(filePath, variableName));
            img_processed = processDataMethod(data, methodName);

            % Ensure frame dimensions are even
            [h, w, ~] = size(img_processed);
            if mod(h, 2) ~= 0
                img_processed = img_processed(1:end-1, :, :);
            end
            if mod(w, 2) ~= 0
                img_processed = img_processed(:, 1:end-1, :);
            end

            % Convert the image data from double [0,1] to uint8 [0,255]
            img_uint8 = uint8(img_processed * 255);

            % Write the frame to the video
            writeVideo(v, img_uint8);

        catch ME
            fprintf('Error processing %s: %s\n', fileName, ME.message);
            continue;
        end
    end
    
    close(v);
    fprintf('Video saved to %s\n', outputVideoFile);
end

function img_processed = processDataMethod(data, methodName)
    switch methodName
        case 'raw_normalized'
            % Raw data normalized
            img_processed = normalizeData(data);
            img_processed = 1 - img_processed; % Invert grayscale
            img_processed = repmat(img_processed', [1, 1, 3]);

        case 'truncated'
            % Truncated IR data between 280 and 292.5 K
            lower_bound = 280; 
            upper_bound = 292.5;
            img_p = data;
            img_p(img_p < lower_bound) = lower_bound;
            img_p(img_p > upper_bound) = upper_bound;
            img_p = (img_p - lower_bound) / (upper_bound - lower_bound);
            img_p = 1 - img_p;
            img_processed = repmat(img_p', [1, 1, 3]);

        case 'highpass_20'
            % High-pass filtering with filterWidth=100
            img_p = applyHighPass(data, 20, false);
            img_processed = repmat(img_p', [1, 1, 3]);

        case 'highpass_100'
            % High-pass filtering with filterWidth=100
            img_p = applyHighPass(data, 100, false);
            img_processed = repmat(img_p', [1, 1, 3]);

        case 'highpass_100_sqrt'
            % High-pass filtering with filterWidth=50 and sqrt enhancement
            img_p = applyHighPass(data, 100, true);
            img_processed = repmat(img_p', [1, 1, 3]);

        case 'highpass_50'
            % High-pass filtering with filterWidth=50 and sqrt enhancement
            img_p = applyHighPass(data, 50, false);
            img_processed = repmat(img_p', [1, 1, 3]);

        case 'highpass_50_sqrt'
            % High-pass filtering with filterWidth=50 and sqrt enhancement
            img_p = applyHighPass(data, 50, true);
            img_processed = repmat(img_p', [1, 1, 3]);

        case 'highpass_20_sqrt'
            % High-pass filtering with filterWidth=20 and sqrt enhancement
            img_p = applyHighPass(data, 20, true);
            img_processed = repmat(img_p', [1, 1, 3]);

        case 'raw_masked'
            % Raw data but mask values above 295K
            maskLimit = 295;
            data(data > maskLimit) = NaN;
            img_p = normalizeDataNaN(data);
            img_p = 1 - img_p; % Invert if desired
            img_processed = repmat(img_p', [1, 1, 3]);

        otherwise
            error('Unknown method: %s', methodName);
    end
end

function img_processed = normalizeData(data)
    % Normalize data to [0,1] based on global min and max
    dmin = min(data(:));
    dmax = max(data(:));
    img_processed = (data - dmin) / (dmax - dmin);
    img_processed(isnan(img_processed)) = 0; % replace NaNs if any
end

function img_processed = normalizeDataNaN(data)
    % Normalize data to [0,1] ignoring NaNs
    dmin = nanmin(data(:));
    dmax = nanmax(data(:));
    img_processed = (data - dmin) / (dmax - dmin);
    % Replace NaNs (masked areas) with 0 after normalization
    img_processed(isnan(img_processed)) = 0;
end

function img_processed = applyHighPass(data, filterWidth, doSqrtEnhance)
    % High-pass filter the data
    lowPass = imgaussfilt(data, filterWidth);
    highPass = data - lowPass;

    % Optionally steepen contrast around zero using sqrt filtering
    if doSqrtEnhance
        highPass = sqrt(abs(highPass)) .* sign(highPass);
    end

    % Normalize using percentiles to avoid outliers
    minVal = prctile(highPass(:), 1);
    maxVal = prctile(highPass(:), 99);
    img_processed = highPass;
    img_processed(img_processed < minVal) = minVal;
    img_processed(img_processed > maxVal) = maxVal;
    img_processed = (img_processed - minVal) / (maxVal - minVal);

    % Invert grayscale (color inversion)
    img_processed = 1 - img_processed;
end
