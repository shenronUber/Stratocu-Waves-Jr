% Import video 
%v = VideoReader("DATA/Baja_202307051401-202307052331_g18_conus_band2_vonkarmans-waves_nolabels.mp4");
%red = squeeze( video(:,:,1,:) );
%size(red)  % 1080, 1920, 115
% frame50=red(:,:,50); % for vonkarmans case

% v = VideoReader("DATA/syn_scwave_warp_modx5.mp4");
%v = VideoReader("DATA/closedcellday_2022_09_06.mp4");
%video = read(v);
%red = squeeze( video(:,:,1,:) ); % red channel brightness only 
%frame =red(100:700,150:880,3); % for synthetic case
%frame2=red(100:700,150:880,4); DT = 30*60; % 4th image in video, 1/2 h later

v = VideoReader("DATA/closedcellday_2022_09_06.mp4");
video = read(v);
red = squeeze( video(:,:,1,:) ); % red channel brightness only 

% Full-res Frames 3 and 4. Flip arrays, images are upside down
fframe = red(1:1000,600:1450,3) ; % for real data image, 3rd image in video
fframe2= red(1:1000,600:1450,4) ; DT = 30*60; % 4th image in video, 1/2 h later
% coordinate arrays at full resolution
fy = 1:size(fframe,1);
fx = 1:size(fframe,2);

% resize images (arrays) to make computations faster  
% shrinkfactor=1;
shrinkfactor=5; invshrinkfactor = 0.2;
frame = imresize(fframe , invshrinkfactor);
frame2= imresize(fframe2, invshrinkfactor);
% coordinate arrays at low resolution
y = (1:size(frame,1) ) *shrinkfactor;
x = (1:size(frame,2) ) *shrinkfactor;


% QUICK SHOW 
figure(1)
image(fx,fy,fframe); colorbar; axis on; hold on; 
contour(x,y,frame2-frame);
title('Frame 1, and contours of Frame2-Frame1')


%%% Wavelet transform inputs (scales and angles) 
% 24 Angles from 0 to pi 
NSCALES = 24;
Angles = 0:pi/NSCALES:pi ;

% a LOGARITHMIC set of 10 Scales, shrinkfactor resizes them
Scales = 10.^(1:.05:1.9);
Scales = Scales/shrinkfactor; 


%%% Call wavelet spectrum 

%%% Ready to test with Morlet 
% wavelet_name = "morlet"; params = [2.0, 3.0, 1.0]
% wavelet_struct = {'name': wavelet_name, 'param': params}
%spec = cwtft2(frame,wavelet="wavelet_name",wavelet_struct,
% scales=Scales,angles=Angles);

% Cauchy one-liner: at two times 
cwt  = cwtft2(frame ,wavelet="cauchy",scales=Scales, angles=Angles);
spec = squeeze( cwt.cfs );

% Second image frame2 
cwt2 = cwtft2(frame2,wavelet="cauchy",scales=Scales, angles=Angles);
spec2= squeeze( cwt2.cfs );

power  = abs(spec) .^2;
power2 = abs(spec2).^2;
xspec  = spec2 .* conj(spec);



% The cross spectrum xspec must be averaged over some samples for its 
% P and Q parts (recast as coherence and angle) to be meaningful 

% Let's hunt for easter eggs of wavelet power, without knowing location. 
% Compute area-averaged power and xspec by angle and scale, away from edges
% where the periodic unpadded unwindowed FFT has shocks (pad/vignette later)

% mean() averages over the first dimension, so two of those will make
% for synthetic data the inner area is smaller
% innerpower = squeeze(  mean(mean( power(100:500,100:600, :,:) ))  );

% for actual closedcell movie data some "inner" box 
% specify it according to array size, max(Scales) away from edges 

buffer = round( max(Scales) );

innerpower = squeeze( mean(mean( power(buffer:size(power,1)-buffer, ...
                                       buffer:size(power,2)-buffer, :,:) )));


% There's a mean increase of power with scale, divide by it 
figure(2)
meanbyscale = squeeze( mean(transpose(innerpower)) );
plot(Scales, meanbyscale); title('mean power by scale'); xlabel('scale (pixels)')

% Normalize by that mean increase with scale, call it anglespec:
anglespec = innerpower .* 0;     % right sized container for angle spectrum
for isc = 1:NSCALES   %size( transpose(Scales) )
    anglespec(:,isc) = squeeze(innerpower(:,isc)) ./ transpose(meanbyscale);
end 

% The angle spectrum, normalized by scale: index coordinates, and actual
figure(3)
pcolor(anglespec); colorbar(); 
xlabel('Angle index'); ylabel('Scale index')
title('areameanpower/meanbyscale')

% This should be polar, but at least it will be stretched like scale & 180
% deg for angle axis label 
figure(4)
pcolor(Angles*180/pi, Scales, anglespec); colorbar(); 
xlabel('Angle (deg)'); ylabel('Scale (pixels, roughly)')
title('areameanpower/meanbyscale (learn to do polar plot)')


% Grab the peaks where anglespec exceeds a threshold strength, and
% illustrate them with a call to image_with_wavelet_overlay
[row,col] = find( imregionalmax(anglespec) & anglespec>1 )
for ipeak = 1:size(row)
    figure(10+ipeak);
    isc = row(ipeak); ian = col(ipeak);
    contours = 2; 
    image_with_wavelet_overlay(frame, spec, Scales, isc,ian, contours); 
    title('scale, angle: '+string( Scales(isc))+' pixels, ' ...
          +string(Angles(ian)*180./pi)+'deg');
end


%%
% The coherence-squared and angle spectra 
innerxspec = squeeze( mean(mean( xspec(buffer:size(xspec,1)-buffer, ...
                                       buffer:size(xspec,2)-buffer, :,:) )));
xpowr = abs( xspec.*conj(spec2) );
inner_coh2 = squeeze( mean(mean( xpowr(buffer:size(xspec,1)-buffer, ...
                                       buffer:size(xspec,2)-buffer, :,:) )))/innerpower;

% Average the complex cross-coherence over spatial area, and extract
% coh,angle from the averaged complex array 
inner_angl = angle(innerxspec);


figure(5)
pcolor(inner_coh2); colorbar(); 
xlabel('Angle index'); ylabel('Scale index')
title('Squared coherence between images')

figure(6)
pcolor(inner_angl); colorbar(); 
xlabel('Angle index'); ylabel('Scale index')
title('Phase angle of xspec averaged over Inner area')


% Peak finder of anglespec should go here, and image with annotation is a
% visual check that the peaks are meaningful as seen by eye 

%figure(7)
%image_with_wavelet_overlay(frame, spec, Scales, 9,8); title('scale 9 angle 8')


function output = image_with_wavelet_overlay(img,spec, Scales, scale, angle,clevfactor)
    % Overlay wavelet power on image 
    image(img); colormap(gray); colorbar; axis on

    hold on

    posLevels = (1:2:9)   /clevfactor;
    negLevels = (-9:2:-1) /clevfactor;

    % Adjust contour levels by Scale as factor (for real/imag), factor^2 (for power)
    sfactor = Scales(scale)
    
    % Real part is crests and trofs, imag is gradients, abs is a magnitude map 
    contour( real(spec(:,:,scale,angle)), LevelList=posLevels*sfactor, EdgeColor='red' );
    contour( real(spec(:,:,scale,angle)), LevelList=negLevels*sfactor, EdgeColor='blue' );
    
    % "power" is amplitude abs() squared 
    contour( (abs(spec(:,:,scale,angle))*sfactor) .^2, EdgeColor='white' );

    % legend
end