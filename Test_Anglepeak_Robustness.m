% Import video 
%v = VideoReader("DATA/Baja_202307051401-202307052331_g18_conus_band2_vonkarmans-waves_nolabels.mp4");
%red = squeeze( video(:,:,1,:) );
%size(red)  % 1080, 1920, 115
% frame50=red(:,:,50); % for vonkarmans case

% v = VideoReader("DATA/syn_scwave_warp_modx5.mp4");
%v = VideoReader("DATA/closedcellday_2022_09_06.mp4");
%video = read(v);
%red = squeeze( video(:,:,1,:) ); % red channel brightness only 
%frame50=red(100:700,150:880,3); % for synthetic case

v = VideoReader("DATA/closedcellday_2022_09_06.mp4");
video = read(v);
red = squeeze( video(:,:,1,:) ); % red channel brightness only 
frame50 = red(1:1000,600:1450,3); % for real data image

% QUICK SHOW 
figure(1)
imshow(frame50); colorbar;
%%% Wavelet transform 

% a set of 7 Angles from 0 to pi (N, NNW, WNW, W, WSW, SSW, )
Angles = 0:pi/24:pi ;

% a LOGARITHMIC set of 10 Scales
% Scales = [2,5,10,20,40,80,160,320,640,1000] % 2,5,10 too small, 1000 too big
% Better range 
Scales = 10.^(1:.05:1.9) ;

cwtCauchy = cwtft2(frame50,wavelet="cauchy",scales=Scales, angles=Angles);
spec = squeeze( cwtCauchy.cfs );


% Let's hunt for easter eggs of wavelet power, without knowing location. 
% Compute area-averaged power by angle and scale: away from the edges 
% mean() averages over the first dimension, so two of those will make

power = abs(spec) .^2;

% for synthetic data (has axes in image area, must clip them) 
% innerpower = squeeze(  mean(mean( power(100:500,100:600, :,:) ))  );

% for actual data (avoid edge effects --> some "inner" box)
innerpower = squeeze(  mean(mean( power(50:800,200:800, :,:) ))  );


% There's a mean increase of power with scale, normalize it away
figure(2)
meanbyscale = squeeze( mean(transpose(innerpower)) );
plot(Scales, meanbyscale); title('mean power by scale'); xlabel('scale (pixels)')

% Normalize by that mean increase with scale, call it anglespec:
anglespec = innerpower .* 0;     % right sized container for angle spectrum
for isc = 1:24   %size( transpose(Scales) )
    anglespec(:,isc) = squeeze(innerpower(:,isc)) ./ transpose(meanbyscale);
end 

% The angle spectrum 
figure(3)
pcolor(anglespec); colorbar(); 
xlabel('Angle index'); ylabel('Scale index')
title('areameanpower/meanbyscale')
hold off 

figure(4)
pcolor(Angles*180/pi, Scales, anglespec); colorbar(); 
xlabel('Angle (deg)'); ylabel('Scale (pixels, roughly)')
title('areameanpower/meanbyscale')
hold off 


figure(7)
image_with_wavelet_overlay(frame50, spec, Scales, 9,8); title('scale 9 angle 8')
figure(8)
image_with_wavelet_overlay(frame50, spec, Scales, 15,14); title('scale 15 angle 14')

% figure(7)
% image_with_wavelet_overlay(frame50, spec, Scales, 8,7); title('scale 8 angle 7')
% figure(8)
% image_with_wavelet_overlay(frame50, spec, Scales, 8,8); title('scale 8 angle 8')
% figure(9)
% image_with_wavelet_overlay(frame50, spec, Scales, 8,9); title('scale 8 angle 9')
% figure(10)
% image_with_wavelet_overlay(frame50, spec, Scales, 8,10); title('scale 8 angle 10')

%%
% Define scales and angles for the upper half
Angles = linspace(0, pi, size(anglespec, 2));  % Angles from 0 to pi (upper half)
Scales = 10.^(1:.05:1.9);  % Logarithmic scales (adjust based on your data)

% Create meshgrid for polar coordinates
[Theta, R] = meshgrid(Angles, Scales);  % Create grid for angles and scales

% Convert polar coordinates to Cartesian for plotting with pcolor
[X, Y] = pol2cart(Theta, R);

% Plot using Cartesian coordinates with pcolor
figure;
pcolor(X, Y, anglespec);  % Use pcolor to plot the wavelet power data
shading interp;  % Smooth shading

% Customize the plot with labels and title
title('Wavelet Power Polar Heatmap');
colorbar;  % Add colorbar to indicate wavelet power

% Remove Cartesian ticks to focus on polar aspect
set(gca, 'XTick', [], 'YTick', []);  
axis equal;  % Ensure equal scaling

% Add detailed angle lines with ticks at each angle
hold on;
angle_ticks = linspace(0, pi, 12);  % Add more granular angle ticks (every 15 degrees)
for angle_rad = angle_ticks
    plot([0 max(X(:))*cos(angle_rad)], [0 max(Y(:))*sin(angle_rad)], 'k--');
    
    % Add text labels for each angle tick (convert to degrees for labeling)
    angle_deg = rad2deg(angle_rad);
    text(max(X(:))*cos(angle_rad)*1.1, max(Y(:))*sin(angle_rad)*1.1, ...
        [num2str(round(angle_deg)), '°'], 'HorizontalAlignment', 'center');
end

% Add radial rings for each scale (more granular)
for i = 1:length(Scales)
    theta_ring = linspace(0, pi, 100);  % Full semicircle for each ring
    [x_ring, y_ring] = pol2cart(theta_ring, Scales(i));
    plot(x_ring, y_ring, 'k--');  % Plot the radial ring
    
    % Add text labels for each scale at the rightmost edge (near 0°)
    text(max(x_ring)*1.05, 0, num2str(Scales(i), '%.2f'), 'HorizontalAlignment', 'left');
end

hold off;

%%

% Define scales and angles for the upper half (0 to pi)
Angles_pos = linspace(0, pi, size(anglespec, 2));  % Angles from 0 to pi (upper half)
Scales = 10.^(1:.05:1.9);  % Logarithmic scales

% Extend angles to cover the full circle (0 to 2*pi)
Angles_full = [Angles_pos, Angles_pos + pi];  % Full angles from 0 to 2*pi

% Create meshgrid for polar coordinates
[Theta, R] = meshgrid(Angles_full, Scales);  % Full angle grid for both halves

% Extend the data to cover the lower half, not just mirrored
anglespec_full = [anglespec, anglespec];  % Repeat the data for the second half (adjust based on physical symmetry)

% Convert polar coordinates to Cartesian for plotting with pcolor
[X, Y] = pol2cart(Theta, R);

% Plot using Cartesian coordinates with pcolor
figure;
pcolor(X, Y, anglespec_full);  % Use pcolor to plot the wavelet power data
shading interp;  % Smooth shading

% Customize the plot with labels and title
title('Wavelet Power Full Rose Plot');
colorbar;  % Add colorbar to indicate wavelet power

% Remove Cartesian ticks 
set(gca, 'XTick', [], 'YTick', []);  
axis equal;  % Ensure equal scaling

% Add angular lines with fractions of pi as labels
hold on;
angle_ticks_pos = linspace(0, pi, 7);  % Upper half ticks (0 to pi)
angle_ticks_neg = linspace(pi, 2*pi, 7); % Lower half ticks (pi to 2*pi)
angle_labels_pos = {'0', '\pi/6', '\pi/3', '\pi/2', '2\pi/3', '5\pi/6', '\pi'};  % Labels for upper half
angle_labels_neg = {'\pi', '7\pi/6', '4\pi/3', '3\pi/2', '5\pi/3', '11\pi/6', '2\pi'};  % Labels for lower half

% Display angle ticks for the upper half
for i = 1:length(angle_ticks_pos)
    angle_rad = angle_ticks_pos(i);
    
    % Upper half
    plot([0 max(X(:))*cos(angle_rad)], [0 max(Y(:))*sin(angle_rad)], 'k--');
    text(max(X(:))*cos(angle_rad)*1.1, max(Y(:))*sin(angle_rad)*1.1, ...
        angle_labels_pos{i}, 'HorizontalAlignment', 'center');
end

% Display angle ticks for the lower half
for i = 1:length(angle_ticks_neg)-1
    angle_rad = angle_ticks_neg(i);
    
    % Lower half
    plot([0 max(X(:))*cos(angle_rad)], [0 max(Y(:))*sin(angle_rad)], 'k--');
    text(max(X(:))*cos(angle_rad)*1.1, max(Y(:))*sin(angle_rad)*1.1, ...
        angle_labels_neg{i}, 'HorizontalAlignment', 'center');
end

% Add radial rings using actual Scales values
ring_radii = Scales;  % Use actual scales for the rings
for i = 1:length(ring_radii)
    if ~ismember(i, [2:3 ,5:10,12:13,15,17])
        theta_ring = linspace(0, 2*pi, 100);  % Full circle for each ring
        [x_ring, y_ring] = pol2cart(theta_ring, ring_radii(i));
        plot(x_ring, y_ring, 'k--');  % Plot the radial ring

        text(max(x_ring)*1.05, 0, num2str(ring_radii(i), '%.2f'), 'HorizontalAlignment', 'left')
    end
end
hold off;


function image_with_wavelet_overlay(img,spec, Scales, scale, angle)
    % Overlay wavelet power on image 
    imshow(img); colorbar; axis on
    
    hold on

    posLevels = 1:2:9;
    negLevels = -9:2:-1;

    % Adjust contour levels by Scale as factor (for real/imag), factor^2 (for power)
    factor = Scales(scale);
    
    % Real part is crests and trofs, imag is gradients, abs is a magnitude map 
    contour( real(spec(:,:,scale,angle)), LevelList=posLevels*factor, EdgeColor='red' );
    contour( real(spec(:,:,scale,angle)), LevelList=negLevels*factor, EdgeColor='blue' );
    
    % "power" is amplitude abs() squared 
    contour( (abs(spec(:,:,scale,angle))*factor) .^2 );

    % legend
end