% Define wavelet transform parameters
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

% Choose a specific scale and angle for plotting
% For instance, select the scale closest to 150 km and an angle of 45 degrees
[~, scale_idx] = min(abs(Scales_km - 50)); % Find scale closest to 20 km
[~, angle_idx] = min(abs(Angles - deg2rad(0))); % Find angle closest to 45 degrees

% Create an impulse image
img_size = 512; % Image resolution (adjust as needed)
impulse_img = zeros(img_size, img_size);
center = floor(img_size / 2) + 1;
impulse_img(center, center) = 1;

% Perform 2D continuous wavelet transform using the 'cauchy' wavelet
cwt_struct = cwtft2(impulse_img, 'wavelet', 'cauchy', 'scales', Scales, 'angles', Angles);

% Extract the wavelet coefficients at the specific scale and angle
coeffs = squeeze(cwt_struct.cfs(:, :, :, scale_idx, angle_idx));

% Separate real and imaginary parts
coeffs_real = real(coeffs);
coeffs_imag = imag(coeffs);

% Generate grid for plotting in km
[x, y] = meshgrid(1:img_size, 1:img_size);
x_shifted = (x - center) * pixel_size_km;
y_shifted = (y - center) * pixel_size_km;

% Identify the distance between two crests empirically
% Extract a horizontal cross-section of the wavelet at the center
cross_section = coeffs_real(center, :);  % Horizontal line at the center
[~, crest_indices] = findpeaks(cross_section);  % Find peaks (crests)

% Calculate crest-to-crest distances in pixels
crest_distances_px = diff(crest_indices);

% Convert crest distances from pixels to km
crest_distances_km = crest_distances_px * pixel_size_km;

% Average crest-to-crest distance (empirical wavelength)
empirical_wavelength_km = mean(crest_distances_km);

% Display results
fprintf('Theoretical Wavelength (km): %.2f\n', Scales_km(scale_idx) * pi / sqrt(2));
fprintf('Empirical Wavelength (km): %.2f\n', empirical_wavelength_km);

% Plot the real part of the wavelet function
figure;
surf(x_shifted, y_shifted, coeffs_real, 'EdgeColor', 'none');
colormap jet;
colorbar;
title(sprintf('Real Part of Wavelet Coefficients\nScale: %.1f km, Angle: %.1f°', Scales_km(scale_idx), rad2deg(Angles(angle_idx))));
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Amplitude');
view(3);

% Plot the imaginary part of the wavelet function
figure;
surf(x_shifted, y_shifted, coeffs_imag, 'EdgeColor', 'none');
colormap jet;
colorbar;
title(sprintf('Imaginary Part of Wavelet Coefficients\nScale: %.1f km, Angle: %.1f°', Scales_km(scale_idx), rad2deg(Angles(angle_idx))));
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Amplitude');
view(3);

% Compute modulus and phase
coeffs_modulus = abs(coeffs);
coeffs_phase = angle(coeffs);

% Plot modulus
figure;
surf(x_shifted, y_shifted, coeffs_modulus, 'EdgeColor', 'none');
colormap hot;
colorbar;
title(sprintf('Modulus of Wavelet Coefficients\nScale: %.1f km, Angle: %.1f°', Scales_km(scale_idx), rad2deg(Angles(angle_idx))));
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Amplitude');
view(3);

% Plot phase
figure;
surf(x_shifted, y_shifted, coeffs_phase, 'EdgeColor', 'none');
colormap hsv;
colorbar;
title(sprintf('Phase of Wavelet Coefficients\nScale: %.1f km, Angle: %.1f°', Scales_km(scale_idx), rad2deg(Angles(angle_idx))));
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Phase (radians)');
view(3);

figure;
plot(x_shifted(center, :), cross_section);
hold on;
plot(x_shifted(center, crest_indices), cross_section(crest_indices), 'ro');
title('Section Transversale avec Crêtes Détectées');
xlabel('Position (km)');
ylabel('Amplitude');
hold off;
