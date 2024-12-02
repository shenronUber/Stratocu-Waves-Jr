%% NEW PROGRAM : Parameters and Setup

% Start timing the entire process
total_time = tic;

% Set path and load images
path = 'C:\Users\admin\Documents\GitHub\Stratocu-Waves-Jr\DATA\2022_09_06\';
image_files = dir(fullfile(path, '*9_06*png'));
images = cell(length(image_files), 1);

% Read and process images
for i = 1:length(image_files)
    fname = fullfile(path, image_files(i).name);
    img = imread(fname);
    red = flipud(img(:, 400:end, 1));
    images{i} = red;
end

% Get image dimensions
[height, width] = size(images{1});
x = 0:(width - 1);
y = 0:(height - 1);
[X, Y] = meshgrid(x, y);

% Pixel size
DX = 2000; % meters per pixel
pixel_size_km = DX / 1000; % kilometers per pixel

% Coordinates in meters
Xm = (X - mean(X(:))) * DX;
%Ym = (Y - mean(Y(:))) * DX;
Ym = (-Y + mean(Y(:))) * DX;

%%

% Wave parameters (first wave)
cphase = 20;            % m/s
wavelength = 150e3;     % meters
direction = 135;        % degrees
zamplitude = 100;       % meters
PBLdepth = 1000;        % m

% Wavenumbers
%k = (2 * pi / wavelength) * sin(direction * pi / 180);
%l = (2 * pi / wavelength) * cos(direction * pi / 180);
% Corrected wave numbers
k = (2 * pi / wavelength) * cos(direction * pi / 180); % Now uses cosine
l = (2 * pi / wavelength) * sin(direction * pi / 180); % Now uses sine


% Amplitude modulation (first wave)
packet_center_x = -400e3;
packet_center_y = -400e3;
packet_width_x = 400e3;
packet_width_y = 300e3;

Ampwindow = exp(-(((Xm - packet_center_x) / packet_width_x).^2 + ...
                  ((Ym - packet_center_y) / packet_width_y).^2));

% Time parameters
time_steps = length(images);
time_resolution = 1800; % seconds between frames

% Initialize the grid to store the synthetic images
grid = zeros(height, width, time_steps);

% Compute angular frequency
omega = cphase * (2 * pi / wavelength);

% Corrected phase calculation
for it = 1:time_steps
    t = (it - 1) * time_resolution;

    % Corrected phase
    phase = k * Xm + l * Ym - omega * t;

    % Vertical displacement
    dz = zamplitude * sin(phase) .* Ampwindow;

    % Horizontal displacements for the first wave
    dxy = (zamplitude / PBLdepth) * wavelength * sin(phase - pi / 2) / DX;

    % Corrected dx and dy calculations
    dx = dxy .* cos(direction * pi / 180) .* Ampwindow; % uses cosine
    dy = dxy .* sin(direction * pi / 180) .* Ampwindow; % uses sine

    % Total displacements
    total_dx = dx;
    total_dy = dy;

    % Get the current image
    img = double(images{it});

    % Coordinates for interpolation
    XI = X - total_dx;
    YI = Y - total_dy;

    % Handle boundaries
    XI = max(min(XI, width), 1);
    YI = max(min(YI, height), 1);

    % Warp the image
    warped_img = interp2(X, Y, img, XI, YI, 'linear', 0);

    % Modulate brightness
    modulated_img = warped_img .* (1 + dz / PBLdepth * 5);

    % Store the synthetic image
    grid(:, :, it) = modulated_img;
end


% Video generation from synthetic images

% Parameters for video creation
video_filename = 'output_video2.mp4'; % Output video filename
frame_rate = 10;                     % Frames per second

% Get the number of frames
num_frames = size(grid, 3);

% Create a VideoWriter object
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = frame_rate;
open(v);

% Generate video frames
figure; % Create a new figure for the video
for i = 1:num_frames
    % Display the frame
    imagesc(grid(:, :, i));
    colormap(gray); % Set grayscale colormap
    axis off;       % Turn off axis
    axis image;     % Maintain aspect ratio
    
    % Add frame to the video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Close the video file
close(v);
disp(['Video saved as ' video_filename]);
