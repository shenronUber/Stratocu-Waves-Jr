% Wavelet transform parameters
NANGLES = 24;
Angles = linspace(0, pi, NANGLES);  % Angles from 0 to pi in radians

% Define scales in km with logarithmic spacing (10 km to 500 km)
min_scale_km = 10;   % Minimum scale in km
max_scale_km = 500;  % Maximum scale in km
NSCALES = 20;        % Number of scales
Scales_km = logspace(log10(min_scale_km), log10(max_scale_km), NSCALES);

% Convert scales to pixels (assuming pixel size is known)
pixel_size_km = 2; % Example: each pixel is 2 km
Scales = Scales_km / pixel_size_km; % Convert scales from km to pixels

% Create an impulse image
img_size = 512; % Image resolution (adjust as needed)
impulse_img = zeros(img_size, img_size);
center = floor(img_size / 2) + 1;
impulse_img(center, center) = 1;

% Prepare arrays to store results
empirical_wavelengths_km = zeros(1, NSCALES);

% Loop over all scales to compute the empirical wavelength
for scale_idx = 1:NSCALES
    % Perform 2D continuous wavelet transform for the current scale
    cwt_struct = cwtft2(impulse_img, 'wavelet', 'cauchy', 'scales', Scales(scale_idx), 'angles', Angles);

    % Extract the wavelet coefficients at the central angle (0 radians)
    coeffs = squeeze(cwt_struct.cfs(:, :, 1, 1));  % First angle (0 radians)

    % Separate real part for crest detection
    coeffs_real = real(coeffs);

    % Extract a horizontal cross-section of the wavelet at the center
    cross_section = coeffs_real(center, :);  % Horizontal line at the center

    % Detect all peaks (crests)
    [~, all_crest_indices] = findpeaks(cross_section);

    % Select only the central peak and the first oscillation on either side
    if ~isempty(all_crest_indices)
        % Find the central peak index (the closest to the impulse location)
        central_index = center;  % Index at the impulse center
        left_crest_idx = all_crest_indices(find(all_crest_indices < central_index, 1, 'last'));
        right_crest_idx = all_crest_indices(find(all_crest_indices > central_index, 1, 'first'));
        
        % Ensure valid indices exist for left and right crests
        if ~isempty(left_crest_idx) && ~isempty(right_crest_idx)
            % Compute the distance between the central peak and one adjacent peak
            crest_distance_px = abs(right_crest_idx - central_index);

            % Convert the distance to km (this represents the full wavelength)
            empirical_wavelengths_km(scale_idx) = crest_distance_px * pixel_size_km;
        else
            empirical_wavelengths_km(scale_idx) = NaN; % Mark as NaN if crests are missing
        end
    else
        empirical_wavelengths_km(scale_idx) = NaN; % Mark as NaN if no peaks are detected
    end
end

% Remove NaN values to avoid issues in fitting
valid_indices = ~isnan(empirical_wavelengths_km);
valid_scales = Scales_km(valid_indices);
valid_wavelengths = empirical_wavelengths_km(valid_indices);

% Compute the slope using linear regression
coefficients = polyfit(valid_scales, valid_wavelengths, 1);
slope = coefficients(1);

% Display the slope
fprintf('Slope of empirical data (wavelength vs scale): %.4f\n', slope);

% Compute theoretical wavelengths for comparison
theoretical_wavelengths_km = Scales_km * pi / sqrt(2);

% Plot scale vs wavelength
figure;
plot(Scales_km, empirical_wavelengths_km, 'ro-', 'LineWidth', 1.5, 'DisplayName', 'Empirical');
hold on;
plot(Scales_km, theoretical_wavelengths_km, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Theoretical');
xlabel('Scale (km)');
ylabel('Wavelength (km)');
title('Scale vs Wavelength for the Cauchy Wavelet');
legend('Location', 'Best');
grid on;
hold off;

% Select the scale index closest to 50 km
[~, scale_idx_50] = min(abs(Scales_km - 50));  % Scale closest to 50 km

% Perform 2D continuous wavelet transform for the selected scale
cwt_struct_50 = cwtft2(impulse_img, 'wavelet', 'cauchy', 'scales', Scales(scale_idx_50), 'angles', Angles);

% Extract the wavelet coefficients at the central angle (0 radians)
coeffs_50 = squeeze(cwt_struct_50.cfs(:, :, 1, 1));  % First angle (0 radians)

% Separate the real part for crest detection
coeffs_real_50 = real(coeffs_50);

% Extract a horizontal cross-section of the wavelet at the center
cross_section_50 = coeffs_real_50(center, :);  % Horizontal line at the center

% Detect peaks (crests) in the cross-section
[all_peak_values, all_crest_indices] = findpeaks(cross_section_50);

% Focus only on the central peak and the first oscillations on either side
central_index = center;  % The index corresponding to the impulse location
left_crest_idx = all_crest_indices(find(all_crest_indices < central_index, 1, 'last'));
right_crest_idx = all_crest_indices(find(all_crest_indices > central_index, 1, 'first'));

% Combine the three indices
selected_crest_indices = [left_crest_idx, central_index, right_crest_idx];
selected_peak_values = cross_section_50(selected_crest_indices);

% Generate x-axis for plotting in km
x_shifted_50 = ((1:img_size) - center) * pixel_size_km;

% Plot the cross-section with detected crests
figure;
plot(x_shifted_50, cross_section_50, 'b-', 'LineWidth', 1.5);  % Cross-section
hold on;
plot(x_shifted_50(selected_crest_indices), selected_peak_values, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);  % Selected crests
grid on;
title(sprintf('Cross-Section with First Oscillation Peaks (Scale: %.1f km)', Scales_km(scale_idx_50)));
xlabel('Position (km)');
ylabel('Amplitude');
hold off;


