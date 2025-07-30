clear; clc; close all;

%% 0. SETUP THE EXPERIMENT ENVIRONMENT

useCumulativeMask   = true; % If true we build one static mask

% --- Load the base frame for the background texture ---
% Use your existing logic to get the desired file.
rootSepacDir = 'C:\Users\admin\Box\GWaves_2023_10_11-14_SEPAC';
instrument   = 'IR';
startDate    = datetime(2023, 10, 12, 14, 15, 0);
endDate      = datetime(2023, 10, 12, 14, 15, 0); % We only need one frame
dataDir = fullfile(rootSepacDir, upper(instrument), 'Data');

% Use your helper function to get the file list (we'll just use the first one)
[fNames, ~, varName] = getDateRangeFiles(dataDir, startDate, endDate);
if isempty(fNames)
    error('No suitable base frame found. Check your path and date settings.');
end
base_frame_path = fullfile(dataDir, fNames{1});
base_frame = double(ncread(base_frame_path, varName));
fprintf('Loaded base frame from: %s\n', fNames{1});

%% 1. DEFINE PARAMETER RANGES TO TEST

% Wave parameters
wavelengths_to_test = [100e3, 150e3, 200e3]; % meters
directions_to_test  = [210, 235, 260];        % degrees
cphases_to_test     = [10, 15, 20];            % m/s

% Advection (drift) parameters
drift_speeds_to_test = [0, 10, 20];           % m/s
drift_angles_to_test = [45, 90, 135, 180, 225]; % degrees

%% 2. INITIALIZE RESULTS TABLE

results = table();
run_counter = 0;
total_runs = numel(wavelengths_to_test) * numel(directions_to_test) * ...
             numel(cphases_to_test) * numel(drift_speeds_to_test) * ...
             numel(drift_angles_to_test);

%% 3. MAIN EXPERIMENT LOOP

fprintf('Starting uncertainty experiment: %d total runs.\n\n', total_runs);

% Nested loops to test every combination
for wl = wavelengths_to_test
    for direc = directions_to_test
        for cp = cphases_to_test
            for ds = drift_speeds_to_test
                for da = drift_angles_to_test
                    
                    run_counter = run_counter + 1;
                    fprintf('--- Running case %d of %d ---\n', run_counter, total_runs);

                    % Call the analysis function for the current combination
                    [perceived_speed, ground_truth_speed, error_ms] = ...
                        measure_synthetic_wave_error(base_frame, wl, direc, cp, ds, da);
                    
                    % Create a result row for this iteration
                    new_row = table(wl, direc, cp, ds, da, ground_truth_speed, perceived_speed, error_ms, ...
                        'VariableNames', {'Wavelength_m', 'Wave_Direction_deg', 'Wave_Speed_ms', ...
                                          'Drift_Speed_ms', 'Drift_Angle_deg', 'Ground_Truth_Speed_ms', ...
                                          'Perceived_Speed_ms', 'Error_ms'});
                                      
                    % Append the new row to the results table
                    results = [results; new_row];
                end
            end
        end
    end
end

fprintf('\n\nExperiment finished.\n');

%% 4. SAVE RESULTS

% Save the complete table to a .mat file for later analysis
save('uncertainty_experiment_results.mat', 'results');
fprintf('Results saved to uncertainty_experiment_results.mat\n');

% Display the head of the results table
disp('Results preview:');
disp(head(results));

% Automatically run the visualization script
visualize_uncertainty_results(results);


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
