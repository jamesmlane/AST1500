%% focus_dm_pokes.m
% Author - James Lane
% Examine poke area of influence in real time.

%% Set local paths
addpath( genpath('../../../src/matlab/pwfs') )

%% Initiating ALPAO
% getenv('ACEROOT')
setenv( 'ACEROOT', 'C:\Users\Admin\Documents\Alpao\AlpaoCoreEngine')
addpath( fullfile(getenv('ACEROOT'), 'matlab') );
acecsStartup();
userStartup
dm.Reset();

% Monitor the DM
dm.StopMonitoring();
dm.StartMonitoring();

% Get current DM commands
load('../assets/flatDMCommands.mat')
% dm.cmdVector(1:dm.nAct)=curDMCommands;

%% Setup the camera
vid = videoinput('pointgrey', 1);
get(vid)
shutter = 20; % Exposure time in ms
flushdata(vid);% clears all frames from buffer
src = getselectedsource(vid);
vid.FramesPerTrigger = 1;
vid.TriggerRepeat = inf;
src.ShutterMode = 'Manual';
src.ExposureMode = 'Off';
src.GainMode = 'Manual';
src.Gain = 0;
src.FrameRateMode = 'Manual';
src.FrameRate = 30;
src.Shutter = shutter;
triggerconfig(vid,'hardware','fallingEdge','externalTriggerMode0-Source0');

%% Initialize the tip-tilt stage
frequency = 50.0; % Hz
amplitude = 50 / 0.7289; % In microns
offset = 2500; % Offset in micro radians
xscale = 0.7289; % Fractional amplitude correction in X
phaseshift = 1.8837;
[E727,Controller] = StartModulation(frequency, amplitude, xscale, phaseshift);

%% Definitions
pauseTime = 0.5;      % Time to pause between command and slope acquisition
nAverage = 5;         % Number of frames to average for each push-pull
nbPushPull = 5;       % Number of push-pull
pushPullValue = 0.12; % Magnitude of push-pull
loopGain = 0.3;       % Loop gain

% Get slopes to check numbers
start(vid)
[slopeX,slopeY] = PWFSGetSlopesModulation(vid);
stop(vid)

% Setup the interaction matrix array
nSlopeX = size(slopeX,1);
nSlopeY = size(slopeY,1);
nSlopes = nSlopeX + nSlopeY;
interactionMatrix = zeros( dm.nAct, nSlopeX+nSlopeY );

%% Create the interaction matrix
wbfig = waitbar(0,'Creating interaction matrix...');

% Record the current DM commands
curDMCommands = dm.cmdVector;

% Start video
start(vid);

% Loop over the number of actuators
for i=1:dm.nAct

  % Setup the array that will hold the push-pull samples
  sDiff = zeros(nSlopes,1);

  % For each actuator loop over the number of push-pull pairs
  for j=1:nbPushPull

    % Acquire slopes
    dm.cmdVector(i)=pushPullValue;
    pause(pauseTime)
    [sPushX,sPushY] = PWFSGetSlopesModulation(vid);
    dm.cmdVector(i)=-pushPullValue;
    pause(pauseTime)
    [sPullX,sPullY] = PWFSGetSlopesModulation(vid);
    dm.cmdVector(i)=curDMCommands(i);

    % Calculate the difference and append
    sDiffX = sPushX-sPullX;
    sDiffY = sPushY-sPullY;
    sDiff = sDiff + [sDiffX;sDiffY];

  end

  % Take the average of the slope vectors
  sVec = sDiff/(2*nbPushPull)/pushPullValue;

  % Append this column into the interaction matrix
  interactionMatrix(i,:) = sVec;
  disp(i);
  waitbar(i/dm.nAct)
end
close(wbfig)

% Stop video
stop(vid)

%% Plot the interaction matrices

f1 = figure('Name','My Interaction Matrix');
xlabel('Lenslet Slope Measurement')
ylabel('Actuator')
imagesc(interactionMatrix);
colorbar

%% Perform SVD and pseudoinverse

[U,S,V] = svd(interactionMatrix);
eigenValues = diag(S);
iS = diag(1./eigenValues);

% Decide on the number of eigenvalues to threshold
condVec = eigenValues(1)./eigenValues;
condNum = 5; % Max Eigenvalue / Eigenvalue threshold
index = condVec > condNum;

% Plot the eigenvalues and the chosen conditioning parameter
f2 = figure('Name','Eigenvalues')
semilogy(eigenValues/eigenValues(1),'.')
line([1,97],[1./condNum,1./condNum],'Color','red')
xlabel('Eigenmodes')
ylabel('Normalized Eigenvalues')

% Threshold the eigenvalues
nThresholded = sum(index);
%nThresholded = 67;
fprintf('%i / 97 modes removed\n',nThresholded)

% Remake the inverse eigenvalue matrix with thresholding
iSth = diag(1./eigenValues(1:end-nThresholded));
[nS,nC] = size(interactionMatrix);
iSth(nC,nS) = 0;

% Compute the command matrix
commandMatrix = V*iSth*U';

%% Plot the command matrices

f3 = figure('Name','My Command Matrix');
imagesc(commandMatrix);
xlabel('Actuator')
ylabel('Pixel Slope Measurement')
colorbar

%% Manual evaluation of the AO loop - flatten the slopes on the PWFS

% Declarations before evaluating the loop
dmIndex = 1:dm.nAct;
start(vid);
f4 = figure('Name','Correcting Slopes');

% Manual loop block
curVector = dm.cmdVector;
[sX,sY] = PWFSGetSlopesModulation(vid);
sVec = [sX;sY];
commandVec = curVector - loopGain*commandMatrix'*sVec;
dm.cmdVector(dmIndex) = commandVec;
[imageData,ts] = PWFSImageCaptureModulation(vid);
imageData = double(imageData);
imagesc(imageData)

% Reset the DM to flat and stop the camera
dm.Reset();
dm.cmdVector(1:dm.nAct) = curDMCommands;
stop(vid)

%% Setup the monitoring of a DM poke

poke_index = 40; % Which actuator to poke

f5 = figure('Name','Monitoring DM Pokes');

ax1 = subplot(3,2,1); % Show the initial X slope measurement
ax2 = subplot(3,2,2); % Show the initial Y slope measurement
ax3 = subplot(3,2,3); % Show the current X slope measurement
ax4 = subplot(3,2,4); % Show the current Y slope measurement
ax5 = subplot(3,2,5); % Show the difference in X
ax6 = subplot(3,2,6); % Show the difference in Y

% Grab slope maps and populate the plots
start(vid);
dm.cmdVector(poke_index) = 0.5;
pause(0.1);
[I1,I2,I3,I4,Sx_init,Sy_init] = PWFSGetSlopeMapsModulation(vid);
dm.cmdVector(poke_index) = curVector(poke_index);
stop(vid);

subplot(ax1);
imagesc(Sx_init);
colorbar;
caxis([-260,260]);
title(ax1,'Initial Sx');

subplot(ax2);
imagesc(Sy_init);
colorbar;
caxis([-260,260]);
title(ax2,'Initial Sy');

% Loop over the number of iterations and update the image each time
n_iteration = 100;
start(vid);
for i = 1:n_iteration

  % Poke the actuator and measure slopes
  dm.cmdVector(poke_index) = 0.5;
  pause(0.1);
  [I1,I2,I3,I4,Sx,Sy] = PWFSGetSlopeMapsModulation(vid);
  dm.cmdVector(poke_index) = curVector(poke_index);

  subplot(ax3);
  imagesc(Sx);
  colorbar;
  caxis([-260,260]);
  title(ax3,'Current Sx');

  subplot(ax4);
  imagesc(Sy);
  colorbar;
  caxis([-260,260]);
  title(ax4,'Current Sy');

  subplot(ax5);
  imagesc(Sx-Sx_init);
  colorbar;
  caxis([-100,100]);
  title(ax5,'Current-Initial Sx');

  subplot(ax6);
  imagesc(Sy-Sy_init);
  colorbar;
  caxis([-100,100]);
  title(ax6,'Current-Initial Sy');

  pause(0.5)
  disp(i)

end
stop(vid);

%% Shut it down
wfs.StopRtd();
wfs.StopSlopeRtd();
wfs.StopAlignmentRtd();
wfs.StopWavefrontRtd();
loop.StopMonitoring();

dm.Reset();
dm.StopMonitoring();

stop(vid);
delete(vid); %closing the connection
clear vid;
clear src;

EndModulation(E727,Controller);
clear E727
clear Controller
