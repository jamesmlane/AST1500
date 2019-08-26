%% pyramid_aberrations.m

% Display varios pyramid wavefront sensor aberrations for reference

% -----------------------------------------------------------------------------

%% Prepare

clear all % Clears all variables
close all % Clear all figures
figDisp = 'yes'; % display all figures? (will slow down the computation)

% Project-relative paths to OOMAO and helper function
addpath(genpath('../../../src/matlab/functions'))
addpath(genpath('../../../src/matlab/OOMAO'))

%% Atmosphere

% r0 = 0.11; %in m, for a site like MMT @ 550 nm
% 
% % Atmosphere object. In V-band, with Fried parameter of r0, with L0 30 meters.
% % Just one single layer for the test bench
% atm = atmosphere(photometry.V,r0,30,...
%     'altitude',0,...
%     'windSpeed',5,...
%     'windDirection',0);

%% WFS parameters

wfsType = 'py'; %'sh' or 'py'
modeType = 'zonal2'; % Modal basis correction - 'zonal2', 'zern', 'KL', 'DH';
% 'zonal1' should be computed for a SHWFS while 'zonal2' is for a PWFS
mag = 10; %guidestar magnitude
Band = 'R'; % guidestar band
dmType = 'reg'; %'ASM' or 'reg', where reg is a DM with square spacing
coupling = 0.35; % deformable mirror actuator coupling
readOutNoise = 2; % wfs detector readout noise

%% Set up the PWFS

nLenslet = 24;
nL   = nLenslet;    % number of lenslets across the aperture (one dimension)
nPx  = 6;           % number of pixels per lenslet
nRes = nL*nPx;      % resolution on the pupil plane (no of pixels)
modulation = 2;     % units of lambda/d
loopGain = 0.5;     % Closed loop integrator gain: generally larger than SHWFS

% Make the PWFS.
% nLenslet is the equivalent size of a SHWFS
% nRes is the number of pixels across the wavefront sensor aperture
% increase modulation to avoid loss of performance due to small linear range
wfs = pyramid(nLenslet,nRes,'modulation',modulation);

%% General parameters

D    = 6.5;                % telescope primary mirror diameter
d    = D/nL;               % lenslet pitch
samplingFreq = 100;        % WFS sampling time
obstructionRatio = 0.0;    % central obscuration ratio of the secondary
fieldOfViewInArcsec = 90;  % field of view in arcseconds

% Makes a telescope. Should reduce nRes so we don't oversample?
tel = telescope(D,'resolution',nRes,...
    'obstructionRatio',obstructionRatio,...
    'fieldOfViewInArcsec',fieldOfViewInArcsec,...
    'samplingTime',1/samplingFreq);

% redefine parameters to match up
nPx = nRes;

%% Definition of a calibration and science source

% Normalize the total flux of the J band to match the R band for a given guide star magnitude
% JL doesn't worry about this
if strcmp(Band,'R')
    ngs = source('wavelength',photometry.R);
else
    ngs = source('wavelength',photometry.J);
end

% And for later use is a science object in K band is instantiated (mainly
% for the strehl computation
science = source('wavelength',photometry.K);

% Propagation of the calibration source to the WFS through the telescope
ngs = ngs.*tel*wfs;

% Initialize the WFS, set the reference slopes
wfs.INIT

%% Show the WFS

% A new frame read-out and slopes computing:
+wfs;
if strcmp(figDisp,'yes')
    % The WFS camera display:
    figure
    subplot(1,2,1)
    imagesc(wfs.camera)
    % The WFS slopes display:
    subplot(1,2,2)
    slopesDisplay(wfs)
end

% -----------------------------------------------------------------------------

%% Define the influence function based on the DM type

switch dmType
    case 'ASM'

        % USE ASM
        load ../../../data/MMT_DM336_Actuators.mat
        nAct = length(ACT);
        actCentrePos(:,1) = ACT(:,1);
        actCentrePos(:,2) = ACT(:,2);
        dmPitch = 6.5/21;
        px = ACT(:,1)*3.25;
        py = ACT(:,2)*3.25;
        % Poke an actuator and measure the true phase offset. What does this
        % Look like??
        bif = gaussianInfluenceFunction(coupling);
        bif.actuatorCoord = (px/dmPitch +1j*py/dmPitch);
        dm = deformableMirror(nAct,'modes',bif,'resolution',tel.resolution,...
            'validActuator',true(1,nAct));

        if strcmp(figDisp,'yes')

            figure(88);
            scatter(px,py)
            title('DM actuator locations')
            axis tight square
            box on

        end

    case 'reg'
        bif = influenceFunction('monotonic',coupling);
        % nAct = nLenslet + 1;
        nAct = 11;
        dm = deformableMirror(nAct,...
            'modes',bif,...
            'resolution',tel.resolution);
end

%% Interaction matrix: DM/WFS calibration

% Now we've defined all of our components. We want to develop a
% relationship between DM commands and phase changes
ngs=ngs.*tel;
stroke = ngs.wavelength/40; %total strength of actuator calibration 'poke'
% This poke is just for calibration. This is because the response
% is linear so you don't have to poke much.

% The choice of mode is how you do your calibration. Do it in terms of
% polynomials (zernike, disk harmonics, etc..), or just poke actuators
switch modeType
    case {'zern','DH','KL','zonal2'}

        %create vector matrix for each mode, scaled by stroke
        radOrd = 25; %set by trying to match at least as many modes as
        %actuators (i.e. degrees of freedom); 25 radial orders of zernikes
        %slightly oversample the number of actuators on the MMT ASM
        j=1;
        kk=0;
        dhRadOrdVec = [];
        for n=0:radOrd
            for k=0:n
                j=j+1;
                dhRadOrdVec = [dhRadOrdVec kk];
            end
            kk = kk+1;
        end

        modeVec = 2:j-1; %exclude piston
        dhRadOrdVec(1) = [];
        zern = zernike(modeVec,'resolution',nRes,'pupil',logical(tel.pupil));
        zernMatrix = pinv(zern.modes)*full(bif.modes);

        % Divide by 1/radial order
        % apparently this is what JP and LAM people do if zernike modes are
        % pushing the PWFS signals beyond the linear regime. I don't think
        % this is an issue for us because we consider a noiseless scenario
        % and so we can use small amplitude 'pokes' or 'shapes'

        % zern.c = zern.c./(zern.n');
        % for jj = 1:length(zern.c)
        %    zern.modes(:,jj) = zern.modes(:,jj)*zern.c(jj);
        % end
        nThresholded = 0;
        switch modeType
            case 'DH'
                %% Use disk harmonics:
                tmp = getDHpoly(radOrd,logical(tel.pupil));
                zern.modes = tmp;
                tmp(:,1) = [];
                zern.modes = tmp;
                %zern.modes = tmp./dhRadOrdVec;
                zernMatrix = pinv(zern.modes)*full(bif.modes);
                nThresholded = 0;
            case 'KL'
                load('MMT_144Res.mat');
                numModes = size(KLi,3);
                modeVec = 1:numModes;
                zern.modes = reshape(KLi,size(KLi,1)^2,size(KLi,3));
                zernMatrix = pinv(zern.modes)*full(bif.modes);
                nThresholded = 0;
            case 'zonal2'
                zernMatrix = eye(nAct);
                modeVec = 1:nAct;
                nThresholded = 1;
        end

        interactionMatrix = zeros(size(wfs.slopes,1),length(modeVec));
        for jj = 1:length(modeVec)

            dm.coefs = stroke*zernMatrix(jj,:)';
            ngs=ngs.*tel*dm*wfs;
            interactionMatrix(:,jj) = wfs.slopes./stroke;
        end
        commandMatrix = pinv(interactionMatrix);
        dmCalib.M = commandMatrix;

        [U,S,V] = svd(interactionMatrix);
        eigenValues = diag(S);

        iS = diag(1./eigenValues(1:end-nThresholded));
        [nS,nC] = size(interactionMatrix);
        iS(nC,nS) = 0;

        % and then the command matrix is derived.
        commandMatrix = V*iS*U';
        dmCalib.M = commandMatrix;
        dmCalib.nThresholded = nThresholded;

    case 'zonal1'

        ngs=ngs.*tel;
        % ?? Where does calibration come from. Comes from OOMAO. This is
        % an object that will help us make our command matrix. Here it
        % pokes the DM and measures the phase changes
        dmCalib = calibration(dm,wfs,ngs,stroke);
        nThresholded = 1;
        dmCalib.nThresholded = nThresholded;
        % OOMAO is doing the inverse for us here.
        commandMatrix = dmCalib.M;

end

% At this point we've made our command / interaction matrices and
% We've defined all the aspects of our telescope.

% -----------------------------------------------------------------------------

%% Prepare to close the AO loop

% Combine the telescope and atmosphere
tel = tel+atm;
% Show the telescope and atmosphere
if strcmp(figDisp,'yes')
    figure(10)
    imagesc(tel)
end

%% Prepare to close the loop

% Flatten the DM
dm.coefs = 0;

% Propagation throught the atmosphere to the telescope
ngs=ngs.*tel;

% Saving the turbulence aberrated phase
turbPhase = ngs.meanRmPhase;

% Propagation to the WFS
ngs=ngs*dm*wfs;

%% Definition of niteration

startDelay = 30;
exposureTime = 60;
nIteration = exposureTime + startDelay;

%% Display turbulence and residual phase

if strcmp(figDisp,'yes')
    figure(11)
    h = imagesc([turbPhase,ngs.meanRmPhase]);
    axis equal tight
    colorbar
    snapnow
end

% Index times for images
u = (0:nIteration-1).*tel.samplingTime;
% set the wavelength of the atmosphere to agree with the guide star
atm.wavelength = ngs.wavelength;

% Piston removed phase variance
totalTheory = phaseStats.zernikeResidualVariance(1,atm,tel);
atm.wavelength = photometry.V;

%% Calculate phase variance in microns and plot

% Phase variance to micron rms converter
% rmsMicron = @(x) 1e6*sqrt(x).*ngs.wavelength/2/pi;
% if strcmp(figDisp,'yes')
%
%     % Where is the actual plotting call??
%     (12)
%     %plot(u,rmsMicron(total),u([1,end]),rmsMicron(totalTheory)*ones(1,2),u,rmsMicron(residue))
%     grid
%     legend('Full','Full (theory)','Residue','Location','Best')
%     xlabel('Time [s]')
%     ylabel('Wavefront rms [\mum]')
% end

%% THEORETICAL PERFORMANCE ANALYSIS - CASE OF THE SINGLE INTEGRATOR LOOP WITH GAIN

varFit      = 0.23*(d/atm.r0)^(5/3)*(atm.wavelength/science(1).wavelength)^2;
varAlias    = 0.07*(d/atm.r0)^(5/3)*(atm.wavelength/science(1).wavelength)^2;
varTempo    = phaseStats.closedLoopVariance(atm, tel.samplingTime,0*tel.samplingTime,loopGain)*(atm.wavelength/science(1).wavelength)^2;
marechalStrehl_lsq_theoretical = 100*exp(-varFit-varAlias-varTempo);

%% WFS noise

magTmp = mag;
ngs.magnitude = magTmp;
ngsTmp = source('wavelength',photometry.R); % Use temporary ngs object to
                                            % normalize photons to R-band
ngsTmp.magnitude = magTmp;
ngs.nPhoton = ngsTmp.nPhoton;

%%

% It can be useful to know the number of photon per subaperture. To do so,
% let separate the atmosphere from the telescope
tel = tel - atm;

%% SCIENCE CAMERA

cam = imager(); % default: 8 pixels across the PSF FWHM at the imaging wavelegth (=4xNyquist)
science = science.*tel*cam;
cam.referenceFrame = cam.frame;
%cam = imager('nyquistSampling',2,'fieldStopSize',1/4*nPx);
+science; % Take a single frame of the science target



%% LOOP INIT

flush(cam) % clear the camera buffer
cam.frame = cam.frame*0; % Reset camera pixels
cam.clockRate    = 1;
cam.exposureTime = exposureTime;

flush(cam)
set(science,'logging',true)
set(science,'phaseVar',[])

%% set the science path

science = science.*tel*dm*cam;

%%

% re-propagate the NGS,
% I think here you reset all the properties of the system but keep the DM
% calibrated so you can just go.
ngs = ngs.*tel*wfs;

%% display the subaperture intensity
if strcmp(figDisp,'yes')
    figure
    intensityDisplay(wfs)
end

%%

% Now the readout noise in photo-electron per pixel per frame rms is set
wfs.camera.readOutNoise = readOutNoise; %in electrons

%%

% Photon-noise is enabled.
wfs.camera.photonNoise = true;

%%

% A pixel threshold is defined. JL doesn't worry
wfs.framePixelThreshold = wfs.camera.readOutNoise;

figure; imagesc(wfs.camera.frame);

tel = tel + atm;

%% number of instances over which the Strehl is computed

totalVar  = zeros(1,nIteration);
residueVar = zeros(1,nIteration);
dm.coefs = 0;

% -----------------------------------------------------------------------------

%% Run the loop

for kIteration=1:nIteration
    % Propagation throught the atmosphere to the telescope, +tel means that
    % all the layers move of one step based on the sampling time and the
    % wind vectors of the layers
    % Having the + here advances the atmosphere by one step and then
    % updates the NGS from start
    ngs=ngs.*+tel;
    +science;

    % Saving the turbulence aberrated phase
    turbPhase = ngs.meanRmPhase;
    % Variance of the atmospheric wavefront
    totalVar(kIteration) = var(ngs);
    % Propagation to the WFS. Here the WFS is updated so it has current
    % slopes
    ngs=ngs*dm*wfs;

    if kIteration == 1
        turbPhaseInit = turbPhase;
    end

    % Variance of the residual wavefront
    residueVar(kIteration) = var(ngs);
    % Computing the DM residual coefficients. Apply the correction using
    % the command matrix, to get residuals, the updating commands.
    residualDmCoefs = commandMatrix*wfs.slopes;
    % Integrating the DM coefficients

    switch modeType
        case {'zonal1','zonal2'}
            % Negative because correction.
            dm.coefs = dm.coefs - loopGain*residualDmCoefs;
        otherwise
            dm.coefs = dm.coefs - (loopGain*((residualDmCoefs')*zernMatrix))';

    end


    if strcmp(figDisp,'yes')

        % Display of turbulence and residual phase
        set(h,'Cdata',[turbPhase,ngs.meanRmPhase])
        drawnow
    end


end

% First simulate the bench with a SH
% Then simulate the bench with a Pyramid
% Want to show a corrected PSF on the bench. Compare bench Strehl with
% theoretical Strehl.
% Simulate what the optimal band pass is going to be.
% Could you have a band pass which is too wide

%%

% Updating the display
if strcmp(figDisp,'yes')

    figure(12)
    text(0.5,1,['SR (Marechal approx):' num2str(marechalStrehl_lsq_theoretical) '%'])
    text(0.5,0.8,['Empirical (Marechal approx):' num2str(100*exp(-(mean(rmsMicron(residueVar(50:nIteration)))*1e-6/ngs.wavelength*2*pi)^2*(atm.wavelength/science(1).wavelength)^2)) '%'])

    set(0,'CurrentFigure',12)
    hold on
    plot(u,rmsMicron(totalVar),'b--',u,rmsMicron(residueVar),'r--')
    legend('Full','Full (theory)','Residue','Full (noise)','Residue (noise)','Location','Best')
end


%% Science camera

cam.strehl

%% Images generated from theory

pupil = logical(turbPhase);
sampling = 4; %pix/fwhm

%% Perfect psf

phase = zeros(size(turbPhase));
wave = pupil.*exp(1i*phase);
H = pupil.*exp(1i*phase).*wave/sqrt(sum(pupil(:).^2));
n = round((length(H)*sampling/10)*10);
HFft = fft2_centre(H,1,n,n);
psf_perf = abs(HFft).^2;
psf_perf = rescale(psf_perf);

[ix,iy] = find(psf_perf == max(psf_perf(:)));
nb = 32; %box width
psf_perf = psf_perf(ix-nb:ix+nb-1,iy-nb:iy+nb-1);

%uncorrected
phase = turbPhaseInit; %in radians @ WFS wavelength
phase = phase*(ngs.wavelength/(2*pi)); %now in m
phase = phase*((2*pi)/science.wavelength); %now in radians @ science wavelength
H = pupil.*exp(1i*phase).*wave/sqrt(sum(pupil(:).^2));
n = round((length(H)*sampling/10)*10);
HFft = fft2_centre(H,1,n,n);
psf_est = abs(HFft).^2;
psf_est = rescale(psf_est).*1e-9;
psf_est = imnoise(psf_est,'poisson');
psf_before = psf_est(ix-nb:ix+nb-1,iy-nb:iy+nb-1);

%corrected
phase = ngs.meanRmPhase;%in radians @ WFS wavelength
phase = phase*(ngs.wavelength/(2*pi)); %now in m
phase = phase*((2*pi)/science.wavelength); %now in radians @ science wavelength
H = pupil.*exp(1i*phase).*wave/sqrt(sum(pupil(:).^2));
n = round((length(H)*sampling/10)*10);
HFft = fft2_centre(H,1,n,n);
psf_est = abs(HFft).^2;
psf_est = rescale(psf_est).*1e-9;
psf_est = imnoise(psf_est,'poisson');
psf_after = psf_est(ix-nb:ix+nb-1,iy-nb:iy+nb-1);

h1 = figure('position',[500 500 1000 500]);
imagesc(log(([psf_before psf_after])))
title('Before and after correction')
h1 = tightfig2(h1);
