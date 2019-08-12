%% focus_dm_pokes.m
% Author - James Lane
% Try and live-view the PWFS while poking the DM to try and test focus
% issues

%% Set local paths
addpath( genpath('../../../src/matlab/pwfs') )

%% Initiating ALPAO
% getenv('ACEROOT')
setenv( 'ACEROOT', 'C:\Users\Admin\Documents\Alpao\AlpaoCoreEngine')
addpath( fullfile(getenv('ACEROOT'), 'matlab') );
acecsStartup();
userStartup
dm.Reset();

%% Define the ALPAO loop

% Switch on the monitoring windows
wfs.StartRtd();
wfs.StartSlopeRtd();
wfs.StartAlignmentRtd();
wfs.StartWavefrontRtd();
wfs.sCam.dit = 0.06; % Set camera exposure time in ms

% Initiate the loop
loop.StartMonitoring();
loop = aceloopLOOP();
loop.set('sWfs', wfs);
loop.set('sWfc', dm);
loop.Online();

%% Setup and monitor the DM
dm.StopMonitoring();
dm.StartMonitoring();
load('./assets/flatDMCommands.m');
dm.cmdVector(1:dm.nAct) = flatDMCommands;

%% Setup the PWFS camera
vid = videoinput('pointgrey', 1);
get(vid)
shutter = 0.01;
flushdata(vid); % clears all frames from buffer
src = getselectedsource(vid);
vid.FramesPerTrigger = 1;
vid.TriggerRepeat = inf;
triggerconfig(vid,'manual');
start(vid);

%% Initialize the tip-tilt stage
frequency = 50.0; % Hz
amplitude = 50 / 0.7289; % In microns
offset = 2500; % Offset in micro radians
xscale = 0.7289; % Fractional amplitude correction in X
phaseshift = 1.8837;
[E727,Controller] = Start_Modulation(frequency, amplitude, xscale, phaseshift);

%% Definitions
pauseTime = 0.1;      % Time to pause between command and slope acquisition
nAverage = 5;         % Number of frames to average for each push-pull
nbPushPull = 5;       % Number of push-pull
pushPullValue = 0.12; % Magnitude of push-pull
loopGain = 0.3;       % Loop gain

% Get slopes to check numbers
[slopeX,slopeY] = PWFSGetSlopes(vid);

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
    [sPushX,sPushY] = PWFSGetSlopes(vid);
    dm.cmdVector(i)=-pushPullValue;
    pause(pauseTime)
    [sPullX,sPullY] = PWFSGetSlopes(vid);
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

%% Create the loop

dmIndex = 1:dm.nAct;
start(vid);
figure();

nsteps = 5000;
for i=1:nsteps

    curVector = dm.cmdVector;

    % Get slopes
    [sX,sY] = PWFSGetSlopesModulation(vid);

    commandVec = curVector - loopGain*commandMatrix'*sVec;

    dm.cmdVector(dmIndex) = commandVec;

    pause(0.01)

    [imageData,ts] = PWFSImageCaptureModulation(vid);
    imageData = double(imageData);
    imagesc(imageData)

    disp(i)
end

%% Setup the monitoring of a DM poke

poke_index = 40; % Which actuator to poke

f5 = figure('Name','Monitoring DM Pokes');

ax1 = subplot(2,2,1) % Show the initial X slope measurement
ax2 = subplot(2,2,2) % Show the initial Y slope measurement
ax3 = subplot(2,2,3) % Show the current X slope measurement
ax4 = subplot(2,2,4) % Show the current Y slope measurement

% Loop over the number of iterations and update the image each time
n_iteration = 1000;
for i = 1:n_iteration




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

End_Modulation(E727,Controller);
clear E727
clear Controller
