magVec = 6:2:14;
%magVec = 6;
for mmm = 1:length(magVec)
    clear atm bif dm ngs ngsTmp science tel wfs zern
    close all
    
    addpath(genpath('./src/matlab/OOMAO_Masen'))
    %% Choose which WFS mode to use
    wfsModel = 'diff'; % Options: 'diff', 'geom'
    
    atm = atmosphere(photometry.V,0.11,30,...
        'altitude',[0,4,10]*1e3,...
        'fractionnalR0',[0.7,0.25,0.05],...
        'windSpeed',[5,10,20],...
        'windDirection',[0,pi/4,pi]);
    
    %% Display figures? This will slow down the simulation
    figDisp = 'no';
    
    %% Simulation parameters
    wfsType = 'py';
    modeType = 'zonal2'; %'zonal2', 'zern', 'KL', 'DH' (use zonal1 for SHWFS + regular DM)
    mag = magVec(mmm);
    %mag = 6; %guidestar magnitude
    Band = 'J'; % 'R' or 'J'
    dmType = 'ASM'; % 'ASM' or 'reg' (regular cartesian)
    switch wfsType
        case 'sh'
            nLenslet = 24;
            nL   = nLenslet;          % number of lenslets
            nPx  = 10;                % number of pixels per lenslet
            nRes = nL*nPx;            % resolution on the pupil plane (no of pixels)
            quadCell = 'true';
            modulation = 'NA'; %units of lambda/d
            
            wfs = shackHartmann(nL,nRes,0.75);
            
            if strcmp(quadCell,'true')
                wfs.camera.resolution = nLenslet*[2 2];
            end
            % Closed loop integrator gain:
            loopGain = 0.3;
            tic
        case 'py'
            nLenslet = 24;
            nL   = nLenslet;        % number of lenslets
            nPx  = 6;               % number of pixels per lenslet
            nRes = nL*nPx;    % resolution on the pupil plane (no of pixels)
            
            modulation = 2;
            if strcmp(Band,'R')
                modTmp = modulation*(1215/640); %match modulation angle between J and R band, R band should be BIGGER to match J band for a given mod angle
            else
                modTmp = modulation;
            end
            
            wfs = pyramid(nLenslet,nRes,'modulation',modTmp,'obstructionRatio',0.10,'minLightRatio',0.05); % increase modulation to avoid loss of performance due to small linear range
            %wfs.binning = 2;
            %wfs = pyramid(nLenslet,nRes,'modulation',2);
            % Closed loop integrator gain:
            loopGain = 0.5;
            
    end
    %% SHWFS parameters
    
    D    = 6.5;                 % telescope primary mirror diameter
    d    = D/nL;              % lenslet pitch
    samplingFreq = 1000;      % WFS sampling time
    obstructionRatio = 0.10;   % central obscuration ratio
    fieldOfViewInArcsec = 90; %fieldOfViewInArcsec
    
    tel = telescope(D,'resolution',nRes,...
        'obstructionRatio',obstructionRatio,'fieldOfViewInArcsec',fieldOfViewInArcsec,'samplingTime',1/samplingFreq);
    
    %redfine parameters to match up
    nPx = nRes;
    
    %% Definition of a calibration source
    % The source class constructor has parameters/value pairs of optional inputs:
    %
    % * the zenith angle [rd] (default: 0rd)
    % * the azimuth angle [rd] (default: 0rd)
    % * the wavelength [m] (default: V band)
    % * the magnitude
    %
    % In the following, an on-axis natural guide star in V band is defined.
    if strcmp(Band,'R')
        ngs = source('wavelength',photometry.R);
    else
        ngs = source('wavelength',photometry.J);
    end
    
    % and for later use a science object in H band is instantiated
    science = source('wavelength',photometry.K);
    
    %%
    % Propagation of the calibration source to the WFS through the telescope
    ngs = ngs.*tel*wfs;
    
    wfs.INIT
    %%
    % A new frame read-out and slopes computing:
    +wfs;
    %%
    if strcmp(figDisp,'yes')
        % The WFS camera display:
        figure
        subplot(1,2,1)
        imagesc(wfs.camera)
        
        %%
        % The WFS slopes display:
        subplot(1,2,2)
        slopesDisplay(wfs)
    end
    
    % GAIN CALIBRATION
    if strcmp(wfsType,'sh')
        wfs.pointingDirection = zeros(2,1);
        pixelScale = ngs.wavelength/...
            (2*d*wfs.lenslets.nyquistSampling);
        tipStep = pixelScale/2;
        nStep   = floor(nPx/nLenslet/3)*2;
        sx      = zeros(1,nStep+1);
        u       = 0:nStep;
        wfs.camera.frameListener.Enabled = false;
        wfs.slopesListener.Enabled = false;
        warning('off','oomao:shackHartmann:relay')
        for kStep=u
            ngs.zenith = -tipStep*kStep;
            +ngs;
            drawnow
            sx(kStep+1) = median(wfs.slopes(1:end/2));
        end
        warning('on','oomao:shackHartmann:relay')
        Ox_in  = u*tipStep*constants.radian2arcsec;
        Ox_out = sx*ngs.wavelength/d/2*constants.radian2arcsec;
        %figure
        plot(Ox_in,Ox_out)
        grid
        slopesLinCoef = polyfit(Ox_in,Ox_out,1);
        wfs.slopesUnits = 1/slopesLinCoef(1);
        
        ngs.zenith = 0;
        wfs.pointingDirection = [];
    end
    
    %% USE ASM
    
    switch dmType
        case 'ASM'
            load ./data/MMT_DM336_Actuators.mat
            nAct = length(ACT);
            actCentrePos(:,1) = ACT(:,1);
            actCentrePos(:,2) = ACT(:,2);
            dmPitch = 6.5/21;
            px = ACT(:,1)*3.25;
            py = ACT(:,2)*3.25;
            %bif = gaussianInfluenceFunction(0.35,dmPitch);
            bif = gaussianInfluenceFunction(0.35);
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
            bif = influenceFunction('monotonic',35/100);
            nAct = nLenslet + 1;
            %         dm = deformableMirror(nAct,'modes',bif,'resolution',tel.resolution,...
            %             'validActuator',true(1,nAct));
            dm = deformableMirror(nAct,...
                'modes',bif,...
                'resolution',tel.resolution,...
                'validActuator',wfs.validActuator);
    end
    
    %% Interaction matrix: DM/WFS calibration
    ngs=ngs.*tel;
    
    %stroke = ngs.wavelength/40; %in m
    stroke = (1215e-9)/40; %fixed to J-band wavelength (so apples to apples comp with visible)
    
    switch modeType
        case {'zern','DH','KL','zonal2'}
            
            %create vector matrix for each mode, scaled by stroke
            radOrd = 25; %set by trying to match at least as many modes as actuators (i.e. degrees of freedom)
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
            modeVec = 1:j-1;
            
            modeVec = 2:j-1; %exclude piston
            dhRadOrdVec(1) = [];
            zern = zernike(modeVec,'resolution',nRes,'pupil',logical(tel.pupil));
            zernMatrix = pinv(zern.modes)*full(bif.modes);
            %% Divide by 1/radial order
            % zern.c = zern.c./(zern.n');
            % for jj = 1:length(zern.c)
            %    zern.modes(:,jj) = zern.modes(:,jj)*zern.c(jj);
            % end
            %%
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
                    load('KL_basis/KL_basis/MMT_144Res.mat');
                    numModes = size(KLi,3);
                    modeVec = 1:numModes;
                    zern.modes = reshape(KLi,size(KLi,1)^2,size(KLi,3));
                    zernMatrix = pinv(zern.modes)*full(bif.modes);
                    nThresholded = 0;
                case 'zonal2'
                    zernMatrix = eye(nAct);
                    modeVec = 1:nAct;
                    nThresholded = 5;
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
            subplot(1,2,2)
            semilogy(eigenValues,'.')
            xlabel('Eigen modes')
            ylabel('Eigen values')
            
            % the n last eigen values are filtered out
            
            condVec = eigenValues(1)./eigenValues;
            condNum = 50;
            index = condVec > condNum;
            %nThresholded = sum(index);
            %nThresholded = (length(modeVec)+1)-336;
            iS = diag(1./eigenValues(1:end-nThresholded));
            [nS,nC] = size(interactionMatrix);
            iS(nC,nS) = 0;
            
            % and then the command matrix is derived.
            commandMatrix = V*iS*U';
            dmCalib.M = commandMatrix;
            dmCalib.nThresholded = nThresholded;
            %dmCalib.nThresholded = 0;
            
            
        case 'zonal1'
            
            ngs=ngs.*tel;
            dmCalib = calibration(dm,wfs,ngs,stroke);
            nThresholded = 5;
            dmCalib.nThresholded = nThresholded;
            commandMatrix = dmCalib.M;
            
    end
    
    %% The closed loop
    % Combining the atmosphere and the telescope
    tel = tel+atm;
    
    %%
    % Resetting the DM command
    dm.coefs = 0;
    %%
    % Propagation throught the atmosphere to the telescope
    ngs=ngs.*tel;
    %%
    % Saving the turbulence aberrated phase
    turbPhase = ngs.meanRmPhase;
    %%
    % Propagation to the WFS
    ngs=ngs*dm*wfs;
    %% Definition of nIteration
    startDelay = 30;
    exposureTime = 60;
    nIteration = exposureTime + startDelay;
    
    %%
    
    if strcmp(figDisp,'yes')
        figure(10)
        imagesc(tel)
        % Display of turbulence and residual phase
        figure(11)
        h = imagesc([turbPhase,ngs.meanRmPhase]);
        axis equal tight
        colorbar
        snapnow
    end
    
    toc
    %snapnow
    u = (0:nIteration-1).*tel.samplingTime;
    atm.wavelength = ngs.wavelength;
    %%
    % Piston removed phase variance
    totalTheory = phaseStats.zernikeResidualVariance(1,atm,tel);
    atm.wavelength = photometry.V;
    %%
    % Phase variance to micron rms converter
    rmsMicron = @(x) 1e6*sqrt(x).*ngs.wavelength/2/pi;
    if strcmp(figDisp,'yes')
        figure(12)
        %plot(u,rmsMicron(total),u([1,end]),rmsMicron(totalTheory)*ones(1,2),u,rmsMicron(residue))
        grid
        legend('Full','Full (theory)','Residue','Location','Best')
        xlabel('Time [s]')
        ylabel('Wavefront rms [\mum]')
    end
    %% THEORETICAL PERFORMANCE ANALYSIS - CASE OF THE SINGLE INTEGRATOR LOOP WITH GAIN
    varFit      = 0.23*(d/atm.r0)^(5/3)*(atm.wavelength/science(1).wavelength)^2;
    varAlias    = 0.07*(d/atm.r0)^(5/3)*(atm.wavelength/science(1).wavelength)^2;
    varTempo    = phaseStats.closedLoopVariance(atm, tel.samplingTime,0*tel.samplingTime,loopGain)*(atm.wavelength/science(1).wavelength)^2;
    marechalStrehl_lsq_theoretical = 100*exp(-varFit-varAlias-varTempo)
    
    % figure(12)
    % text(0.5,1,['SR (Marechal approx):' num2str(marechalStrehl_lsq_theoretical) '%'])
    % text(0.5,0.8,['Empirical (Marechal approx):' num2str(...
    %     100*exp(-(mean(rmsMicron(residue(50:end)))*1e-6/ngs.wavelength*2*pi)^2*(atm.wavelength/science(1).wavelength)^2)) '%'])
    
    
    
    %% WFS noise
    % Noise can be added to the wavefront sensor but first we need to set the
    % star magnitude.
    magTmp = mag;
    ngs.magnitude = magTmp;
    ngsTmp = source('wavelength',photometry.R); %normalize photons to R-band
    ngsTmp.magnitude = magTmp;
    ngs.nPhoton = ngsTmp.nPhoton;
    
    %%
    % It can be useful to know the number of photon per subaperture. To do so,
    % let separate the atmosphere from the telescope
    tel = tel - atm;
    
    %% SCIENCE CAMERA
    cam = imager();
    science = science.*tel*cam;
    cam.referenceFrame = cam.frame;
    +science;
    
    flush(cam)
    cam.frame = cam.frame*0;
    cam.clockRate = 1;
    cam.exposureTime = exposureTime;
    
    flush(cam)
    set(science,'logging',true)
    set(science,'phaseVar',[])
    
    %% set the science path
    science = science.*tel*dm*cam;
    %%
    % re-propagate the source,
    ngs = ngs.*tel*wfs;
    %%
    % and display the subaperture intensity
    if strcmp(figDisp,'yes')
        figure
        intensityDisplay(wfs)
    end
    %%
    % Now the readout noise in photo-electron per pixel per frame rms is set
    wfs.camera.readOutNoise = 0;
    %%
    % Photon-noise is enabled.
    wfs.camera.photonNoise = true;
    %%
    % A pixel threshold is defined
    wfs.framePixelThreshold = wfs.camera.readOutNoise;
    %%
    nInst = 10;
    strehlVecMaster = zeros(1,nInst);
    for L = 1:nInst
        strehlVec = zeros(1,length(1:1:nIteration));
        %%
        % The loop is closed again
        total  = zeros(1,nIteration);
        residue = zeros(1,nIteration);
        dm.coefs = 0;
        tel = tel + atm;
        for kIteration=1:nIteration
            % Propagation throught the atmosphere to the telescope, +tel means that
            % all the layers move of one step based on the sampling time and the
            % wind vectors of the layers
            ngs=ngs.*+tel;
            % Saving the turbulence aberrated phase
            turbPhase = ngs.meanRmPhase;
            % Variance of the atmospheric wavefront
            total(kIteration) = var(ngs);
            % Propagation to the WFS
            ngs=ngs*dm*wfs;
            % Variance of the residual wavefront
            residue(kIteration) = var(ngs);
            % Computing the DM residual coefficients
            residualDmCoefs = commandMatrix*wfs.slopes;
            % Integrating the DM coefficients
            
            switch modeType
                case {'zonal1','zonal2'}
                    dm.coefs = dm.coefs - loopGain*residualDmCoefs;
                otherwise
                    dm.coefs = dm.coefs - (loopGain*((residualDmCoefs')*zernMatrix))';
                    
            end
            
            strehlVec(1,kIteration) = exp(-(rmsMicron(residue(kIteration))*((2*pi)/(science.photometry.wavelength*1e6))).^2);
            
            if strcmp(figDisp,'yes')
                % Display of turbulence and residual phase
                set(h,'Cdata',[turbPhase,ngs.meanRmPhase])
                drawnow
            end
        end
        %%
        % Updating the display
        if strcmp(figDisp,'yes')
            figure(12)
            text(0.5,1,['SR (Marechal approx):' num2str(marechalStrehl_lsq_theoretical) '%'])
            text(0.5,0.8,['Empirical (Marechal approx):' num2str(...
                100*exp(-(mean(rmsMicron(residue(50:end)))*1e-6/ngs.wavelength*2*pi)^2*(atm.wavelength/science(1).wavelength)^2)) '%'])
            
            set(0,'CurrentFigure',12)
            hold on
            plot(u,rmsMicron(total),'b--',u,rmsMicron(residue),'r--')
            legend('Full','Full (theory)','Residue','Full (noise)','Residue (noise)','Location','Best')
            
            
            figure; plot(strehlVec(1,:),'*')
        end
        Strehl_marechal = mean(strehlVec(1,end-49:end))
        
        strehlVecMaster(1,L) = Strehl_marechal;
        
    end
    
    str = strcat('./',num2str(wfsType),'_',num2str(modeType),'_nThresh',num2str(nThresholded),'_',num2str(magTmp),'mag_',Band,'band','_modulation',num2str(modulation));
    
    save(str,'strehlVecMaster')
end
