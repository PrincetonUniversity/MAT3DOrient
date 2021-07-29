%% This script demonstrates how one may use the package MAT3DOrient
% in this case, we will simulate the spectral response of a 92nmx40nm Au
% nanorod, and then calculate the Cramer-Rao lower bounds of precision in
% theta and phi at two different photon flux levels, then plot these
% results.

%% first add common functions
addpath('CommonFunctions');  % common functions are in folder above

CRs = CramerRaoFunctions; % initialise our Cramer-Rao Functions
NanoGen = NanoparticleFunctions; % intialise our Nanoparticle simulation functions

thetas = linspace(deg2rad(0), deg2rad(90), 181); % theta range
phis = linspace(deg2rad(0), deg2rad(180), 361); % phi range

time = 1; % time is 1 of intensity unit, therefore when we add in intensities we simply are adding photon counts

Itot = ([10, 100, 1000, 10000]); % we will simply explore the effect of changing intensity by three orders of magnitude
n_detectors = 4; % we will assume four detectors

% General nanoparticle parameters
element = 'Au'; % we will think of an Au rod
shape = 'Rod'; % it is a rod
AHFactor = 2.8; % this ad-hoc factor corrects for the fact that the theory of Yu et al has a weaker transverse resonance than found experimentally
wvl = 200:0.01:2000; % set up wavelengths
L = 92; % 92 nm length rod
Width = 40; % 40 nm length rod
R = L./Width; % get R factor
T = 293.15; % temperature
cv = 0; % glycerol fraction

% Experimentally-defined constants
NACond = 1.3; % condenser NA
NAObj = 0.7; % objective NA

[exttotal, extlongitudinal, exttransverse, scattertotal, ~, ~, ~, ~, ~]...
    = NanoGen.YuSpectra(wvl,shape,T,cv,L,R,element,AHFactor); % get spectra
[~, maxind] = max(extlongitudinal); % get central plasmon wavelength
centrewavelength = wvl(maxind); % this will be central wavelength of our ``lamp''


FWHM = 2.*sqrt(2.*log(2)); % FWHM parameter
lampwidth = 1; % lamp width in nm

wavelengths = linspace(centrewavelength-(3*lampwidth./FWHM), centrewavelength+(3*lampwidth./FWHM), 1000); % our lamp is 1000 data points
IW = @(wvl) normpdf(wvl,centrewavelength,lampwidth./FWHM)./max(normpdf(wvl,centrewavelength,lampwidth./FWHM)); % our lamp is a gaussian centred on the central wavelength
    
[a11, a13, a33] = NanoGen.Yuavalues(wavelengths,IW,T,cv,L,R,element,2.8,'Rod'); %ad-hoc 2.8 factor based on comparison with xpt    
n0 = CRs.n_m(wavelengths, T, cv); % get refactive index
[~, ~, A, B, C, H] = CRs.InstrResp(NACond, NAObj, n0); % get factors needed for thetadepf
thetadepf = CRs.thetadepgetscatter(a11, a13, a33, A, B, H); % get thetadepf

CramerRaoMThetas = zeros(length(thetas), length(phis), length(Itot));
CramerRaoMPhis = zeros(length(thetas), length(phis), length(Itot));
for i = 1:length(Itot)
    warning('off','all');
    fprintf('Intensity %d Photons\n',Itot(i))
    I = Itot(i);
    IBeta = I;
    IBeta = repmat(IBeta, n_detectors, 1);
    [CramerRaoMThetas(:,:,i), CramerRaoMPhis(:,:,i)] = CRs.CramerRaoScatterNB(time, thetas', phis', IBeta,...
        wavelengths, T, cv, NAObj, NACond, a11, a13, a33, n_detectors, thetadepf);
end


% plot
figure(1)
for k1 = 1:4
    subplot(2,2,k1)
    contourf(rad2deg(phis),rad2deg(thetas),rad2deg(CramerRaoMThetas(:,:,k1)));
    xlabel(['\phi'])
    ylabel(['\theta'])
end
hp4 = get(subplot(2,2,4),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.02  hp4(2)+hp4(3)*2.1])

% plot
figure(2)
for k1 = 1:4
    subplot(2,2,k1)
    contourf(rad2deg(phis),rad2deg(thetas),rad2deg(CramerRaoMPhis(:,:,k1)));
    xlabel(['\phi'])
    ylabel(['\theta'])
end
hp4 = get(subplot(2,2,4),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.02  hp4(2)+hp4(3)*2.1])