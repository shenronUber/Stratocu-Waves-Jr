%% Parameters and Setup

% Define access token and base URL
accessToken = 'MG9Uch8NWvi2dB7ApfoZkVBGCSDtz9lo';  % Replace with your actual developer token
baseURL = 'https://api.box.com/2.0/';

% Set up web options with headers for authorization
options = weboptions('HeaderFields', {'Authorization', ['Bearer ', accessToken]});

% Define parameters
rootFolderID = '288122029977';  % Replace with your root folder ID
dataFolderName = '4km_SEPAC_VIS';  % Choose '4km_SEPAC_IR' or '4km_SEPAC_VIS'
startDate = datetime(2023, 11, 1);
endDate = datetime(2023, 11, 7);
downloadDir = 'downloaded_data';
outputVideoFile = 'output_video.mp4';

% Retrieve file list
%fileList = retrieveFilesForPeriod(rootFolderID, dataFolderName, startDate, endDate, options, baseURL);

% Download files
%downloadFiles(fileList, downloadDir, options, baseURL);
%%
% Only call the video creation function if files are already downloaded
createVideoFromFiles(downloadDir, 'IR');

%% Functions

function fileList = retrieveFilesForPeriod(rootFolderID, dataType, startDate, endDate, options, baseURL)
    % Initialize file list
    fileList = {};

    % Navigate to the folder for the data type
    dataFolderID = getFolderIDByName(rootFolderID, dataType, options, baseURL);
    if isempty(dataFolderID)
        error('Data type folder "%s" not found in root folder.', dataType);
    end

    % Get contents of the data folder
    folderContents = getFolderContents(dataFolderID, options, baseURL);
    entriesArray = folderContents.entries;  % Array of entry structs

    % Loop through each struct in the entries array
    for j = 1:numel(entriesArray)
        entry = entriesArray{j};  % Extract each struct array in entries cell

        % Now loop through the struct array in this entry
        for i = 1:numel(entry)
            fileEntry = entry(i);  % Access individual file entry struct
            
            if strcmp(fileEntry.type, 'file')
                % Extract the timestamp from the filename to filter by date
                fileName = fileEntry.name;
                
                % Remove the extra .4km.nc extension if present
                cleanedFileName = regexprep(fileName, '\.nc\.4km$', '');
                
                % Use regex to find the 14-character timestamp after "_s"
                timestampStr = regexp(cleanedFileName, '_s(\d{14})_', 'tokens', 'once');
                
                % Ensure the timestamp string is correctly formatted
                if ~isempty(timestampStr)
                    % Parse the timestamp with 14 characters, including extra seconds
                    fileTime = datetime(timestampStr{1}, 'InputFormat', 'uuuuDDDHHmmssSS');
                    
                    % Check if file falls within the start and end dates
                    if fileTime >= startDate && fileTime <= endDate
                        fileInfo = struct('name', fileEntry.name, 'id', fileEntry.id);
                        fileList{end+1} = fileInfo;
                    end
                else
                    warning('Unexpected timestamp format in file name: %s', fileName);
                end

            end
        end
    end
end

function folderID = getFolderIDByName(parentFolderID, folderName, options, baseURL)
    contents = getFolderContents(parentFolderID, options, baseURL);
    entries = contents.entries;

    folderID = '';

    for i = 1:length(entries)
        entry = entries{i};

        if iscell(entry)
            entry = entry{1};
        end

        if isstruct(entry) && isfield(entry, 'type') && isfield(entry, 'name')
            if strcmp(entry.type, 'folder') && strcmp(entry.name, folderName)
                folderID = entry.id;
                return;
            end
        else
            % Display error message with available fields if unexpected structure
            error('Entry does not have expected fields "type" or "name". Available fields are: %s', strjoin(fieldnames(entry), ', '));
        end
    end
end

function folderContents = getFolderContents(folderID, options, baseURL)
    offset = 0;
    limit = 1000;  % Box API limit per request
    allEntries = {};  % Cell array to accumulate all entries

    while true
        % Build the URL with pagination parameters for each request
        paginatedURL = sprintf('%sfolders/%s/items?limit=%d&offset=%d', baseURL, folderID, limit, offset);
        
        % Fetch the current batch of entries
        contentBatch = webread(paginatedURL, options);

        % Append current entries to the allEntries cell array
        allEntries = [allEntries; contentBatch.entries];

        % If the batch is less than the limit, we've reached the end
        if numel(contentBatch.entries) < limit
            break;
        end

        % Otherwise, increment the offset for the next batch
        offset = offset + limit;
    end

    % Combine all entries into folderContents structure
    folderContents.entries = allEntries;
end

function downloadFiles(fileList, downloadDir, options, baseURL)
    if ~exist(downloadDir, 'dir')
        mkdir(downloadDir);
    end

    for i = 1:length(fileList)
        fileInfo = fileList{i};
        fileName = fileInfo.name;
        fileID = fileInfo.id;
        fprintf('Downloading %s...\n', fileName);

        % Construct the download URL for the file
        downloadURL = [baseURL, 'files/', fileID, '/content'];
        
        % Define the path to save the file locally
        filePath = fullfile(downloadDir, fileName);

        % Try downloading the file with websave
        try
            websave(filePath, downloadURL, options);
            fprintf('Successfully downloaded %s\n', fileName);
        catch ME
            fprintf('Failed to download %s: %s\n', fileName, ME.message);
        end
    end
end

function createVideoFromFiles(downloadDir, dataType)
    % dataType should be 'IR' or 'VIS'

    % List all .nc files in the download directory
    ncFiles = dir(fullfile(downloadDir, '*.nc*'));  % Adjusted to include .nc and .nc.4km files

    if isempty(ncFiles)
        error('No .nc files found in the download directory.');
    end

    % Initialize fileTimestamps as an empty datetime array
    fileTimestamps = datetime([], 'ConvertFrom', 'datenum');
    fileNames = {};
    
    for i = 1:length(ncFiles)
        fileName = ncFiles(i).name;
        timestampStr = extractBetween(fileName, '_s', '_e');
        if isempty(timestampStr)
            continue;
        end
    
        % Use the correct timestamp format
        timestampFormat = 'uuuuDDDHHmmssSS';  % Adjusted to match the timestamp format in the filenames
        try
            fileTimestamp = datetime(timestampStr{1}, 'InputFormat', timestampFormat);
            fileTimestamps(end+1) = fileTimestamp;  % This now works as fileTimestamps is a datetime array
            fileNames{end+1} = fileName;
        catch
            continue;
        end
    end

    % Sort files based on extracted timestamps
    [fileTimestamps, sortIdx] = sort(fileTimestamps);
    fileNames = fileNames(sortIdx);

    % Now, process the files depending on dataType
    if strcmpi(dataType, 'IR')
        % For IR data, create one video covering the entire time range
        startTime = datestr(fileTimestamps(1), 'yyyymmddHHMMSS');
        endTime = datestr(fileTimestamps(end), 'yyyymmddHHMMSS');
        outputVideoFile = sprintf('IR_%s_%s.mp4', startTime, endTime);

        % Create video
        createVideo(fileNames, fileTimestamps, downloadDir, outputVideoFile);

    elseif strcmpi(dataType, 'VIS')
        % Shift each timestamp to the beginning of the day to get unique dates
        uniqueDates = unique(dateshift(fileTimestamps, 'start', 'day'));
        for d = 1:length(uniqueDates)
            % Find indices where the timestamp matches the current unique date
            dayFilesIdx = find(dateshift(fileTimestamps, 'start', 'day') == uniqueDates(d));
            dayFileNames = fileNames(dayFilesIdx);
            dayFileTimestamps = fileTimestamps(dayFilesIdx);
            dateStr = datestr(uniqueDates(d), 'yyyymmdd');
            outputVideoFile = sprintf('VIS_%s.mp4', dateStr);
        
            % Create video for this day
            createVideo(dayFileNames, dayFileTimestamps, downloadDir, outputVideoFile);
        end
    else
        error('Invalid dataType. Must be ''IR'' or ''VIS''.');
    end
end

function createVideo(fileNames, fileTimestamps, downloadDir, outputVideoFile)
    % Initialize video writer
    v = VideoWriter(outputVideoFile, 'MPEG-4');
    v.FrameRate = 2;  % Adjust frame rate as needed
    open(v);

    for i = 1:length(fileNames)
        fileName = fileNames{i};
        filePath = fullfile(downloadDir, fileName);
        fprintf('Processing %s...\n', fileName);

        try
            % Check the file for the variable of interest
            info = ncinfo(filePath);
            if any(strcmp({info.Variables.Name}, 'CMI'))
                variableName = 'CMI';
            elseif any(strcmp({info.Variables.Name}, 'Rad'))
                variableName = 'Rad';
            else
                % List all available variables in the file for reference
                availableVars = {info.Variables.Name};
                warning('No suitable variable found in file %s. Available variables: %s. Skipping this file.', ...
                        fileName, strjoin(availableVars, ', '));
                continue;
            end

            % Read the selected data
            data = ncread(filePath, variableName);

            % Display the data without any additional elements
            % Normalize data to the range [0, 1] for image display
            data = double(data);
            data = (data - min(data(:))) / (max(data(:)) - min(data(:)));

            % Create an RGB image from the grayscale data
            img = repmat(data', [1, 1, 3]);  % Transpose data to match orientation

            % Write the frame to the video
            writeVideo(v, img);

        catch ME
            fprintf('Error processing %s: %s\n', fileName, ME.message);
            continue;
        end
    end

    close(v);
    fprintf('Video saved to %s\n', outputVideoFile);
end
