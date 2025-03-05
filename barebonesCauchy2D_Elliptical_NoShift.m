function cwtstruct = barebonesCauchy2D_Elliptical_NoShift(X, scales, angles, ...
    coneAngle, sigmaX, sigmaY, alpha, L, M, normType)
% BAREBONESCauchy2D_Elliptical_NoShift 2D CWT with elliptical Cauchy wavelet, no scale shift.
%
%   cwtstruct = barebonesCauchy2D_Elliptical_NoShift(X, SCALES, ANGLES, ...)
%   returns a transform similar to cwtft2(...,'wavelet','cauchy'), but:
%     - k0=0 (no frequency shift)
%     - elliptical envelope in frequency
%     - alpha controls the "radial" decay (bandwidth)
%
%   Inputs:
%     X         : 2D real matrix (e.g. grayscale image).
%     SCALES    : Vector of scales, must be >= 1 (e.g. 2.^(0:5)).
%     ANGLES    : Vector of angles in radians (e.g. 0:pi/8:pi).
%     coneAngle : Cauchy wavelet parameter (default: pi/6).
%     sigmaX    : Elliptical coefficient along ωX (default: 1).
%     sigmaY    : Elliptical coefficient along ωY (default: 1).
%     alpha     : Overall radial decay factor (NEW parameter).
%     L, M      : Polynomial exponents (default: 4,4).
%     normType  : 'L0', 'L1', or 'L2' (default: 'L2').
%
%   Output: a struct with fields:
%       - cfs      : wavelet coefficients, size (H, W, #scales, #angles)
%       - scales   : the input scales
%       - angles   : the input angles
%       - wav      : struct with .wname='cauchy', .param={coneAngle, sigmaX, sigmaY, alpha, L, M}
%       - wav_norm : normalization factors
%       - meanSIG  : mean of X
%
%   Example usage:
%       X = phantom(256);
%       scales = 2.^(0:5);
%       angles = 0:pi/6:pi;
%       coneAngle = pi/6; 
%       sigmaX = 0.5;  sigmaY = 1.5; 
%       alpha  = 2;   % increase alpha => narrower radial passband
%       L=4; M=4;
%       cwtstruct = barebonesCauchy2D_Elliptical_NoShift( ...
%           X, scales, angles, coneAngle, sigmaX, sigmaY, alpha, L, M );
%       % Compare results to standard cwtft2(...) for reference

%--------------------------------------------------------------------------
% 1) Handle defaults
if nargin < 10 || isempty(normType), normType = 'L2'; end
if nargin < 9 || isempty(M),         M = 4;           end
if nargin < 8 || isempty(L),         L = 4;           end
if nargin < 7 || isempty(alpha),     alpha = 1;       end
if nargin < 6 || isempty(sigmaY),    sigmaY = 1;      end
if nargin < 5 || isempty(sigmaX),    sigmaX = 1;      end
if nargin < 4 || isempty(coneAngle), coneAngle = pi/6;end
if nargin < 3
    error('At minimum: need X, scales, angles');
end

validateattributes(scales,{'double','single'},{'vector','>=',1,'nonempty'});
validateattributes(angles,{'double','single'},{'vector','nonempty'});

switch upper(normType)
    case 'L2', normPOW = 1;
    case 'L1', normPOW = 0;
    case 'L0', normPOW = 2;
    otherwise
        error('normType must be ''L0'', ''L1'', or ''L2''.');
end
pipow = 0.5/(2*pi);  

%--------------------------------------------------------------------------
% 2) Basic info & FFT
X = double(X);
[H,W] = size(X);
meanSIG = mean(X(:));
fimg = fft2(X);

%--------------------------------------------------------------------------
% 3) Build frequency grid
W2  = floor((W-1)/2);
H2  = floor((H-1)/2);
W_puls = (2*pi/W) * [0:W2, (W2 - W + 1):-1];
H_puls = (2*pi/H) * [0:H2, (H2 - H + 1):-1];
[xx, yy] = meshgrid(W_puls, H_puls);
dxx = xx(1,2) - xx(1,1);
dyy = yy(2,1) - yy(1,1);
dxxdyy = abs(dxx*dyy);

nbSca = numel(scales);
nbAng = numel(angles);

cfs = zeros(H, W, nbSca, nbAng);
wav_norm = zeros(nbSca, nbAng);

%--------------------------------------------------------------------------
% 4) Main scale-angle loop
for sIdx = 1:nbSca
    valSca = scales(sIdx);
    for aIdx = 1:nbAng
        valAng = angles(aIdx);
        % Rotate + scale freq grid
        nxx = valSca * ( cos(valAng)*xx - sin(valAng)*yy );
        nyy = valSca * ( sin(valAng)*xx + cos(valAng)*yy );
        
        % Evaluate elliptical wavelet w/ no shift, extra alpha
        wft = cauchyWavelet2DElliptical_NoShift(nxx, nyy, ...
            coneAngle, sigmaX, sigmaY, alpha, L, M);

        % Multiply by scale^normPOW for normalization
        mask = valSca^normPOW * wft;

        % ifft2 => cfs
        cfs(:,:,sIdx,aIdx) = ifft2( fimg .* conj(mask) );

        % wavelet normalization factor
        maskEnergy = sum(abs(mask(:)).^2);
        wav_norm(sIdx,aIdx) = (maskEnergy * dxxdyy)^pipow;
    end
end

%--------------------------------------------------------------------------
% 5) Build output struct
WAV = struct('wname','cauchy',...
             'param',{{coneAngle, sigmaX, sigmaY, alpha, L, M}});
cwtstruct = struct('cfs',     cfs,...
                   'scales',  scales,...
                   'angles',  angles,...
                   'wav',     WAV,...
                   'wav_norm',wav_norm,...
                   'meanSIG', meanSIG);
end

%====================================================================
function wft = cauchyWavelet2DElliptical_NoShift(omegaX, omegaY, ...
    coneAngle, sigmaX, sigmaY, alpha, L, M)
% Elliptical Cauchy wavelet in freq domain, no shift, plus alpha for radial control

% 1) Polynomial factor
dot1  = sin(coneAngle)*omegaX + cos(coneAngle)*omegaY;
dot2  = -sin(coneAngle)*omegaX + cos(coneAngle)*omegaY;
coeff = (dot1.^L) .* (dot2.^M);

% 2) No shift => k0=0
k0 = 0;

% 3) Elliptical + radial factor alpha
%    rad2 = alpha * [0.5*(sigmaX*(omegaX-k0)^2 + sigmaY*(omegaY)^2 ) ]
%    => multiply everything by alpha so bigger alpha => narrower passband
rad2 = 0.5 * alpha .* ( sigmaX*(omegaX - k0).^2 + sigmaY*(omegaY).^2 );

% 4) "pond" mask
pond = tan(coneAngle)*omegaX > abs(omegaY);

% 5) Combine
wft = pond .* coeff .* exp(-rad2);

% 6) Normalize so max(abs(wft)) == 1
mx = max(abs(wft(:)));
if mx > 0
    wft = wft / mx;
end
end
