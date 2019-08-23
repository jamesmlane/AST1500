%% run_pwfs.m
% Author - James Lane
% Run the Pyramid Wavefront Sensor with modulation

%% Set local paths
addpath( genpath('../../src/matlab/pwfs') )

%% Setup the camera
imaqreset
vid = videoinput('pointgrey', 1); 
get(vid)
shutter = 10; % Exposure time in ms
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
load('./assets/flatDMCommands.mat')
% dm.cmdVector(1:dm.nAct)=curDMCommands;

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

%% Initialize the tip-tilt stage
frequency = 100.0; % Hz
xscale = 0.795; % Fractional amplitude correction in X
phaseshift = 1;
amplitude = 450 / xscale; % In microns
offset = 2500; % Offset in micro radians

[E727,Controller] = StartModulation(frequency, amplitude, xscale, phaseshift);

%% Take background images
start(vid)
pause(0.1)
background = TakeBackgroundImage(vid);
stop(vid)

%% Definitions
pauseTime = 0.1;      % Time to pause between command and slope acquisition
nAverage = 5;         % Number of frames to average for each push-pull
nbPushPull = 5;       % Number of push-pull
pushPullValue = 0.12; % Magnitude of push-pull
loopGain = 0.1;       % Loop gain

% Get slopes to check numbers
start(vid)
[slopeX,slopeY] = PWFSGetSlopesModulation(vid,background);
stop(vid)

% Setup the interaction matrix array
nSlopeX = size(slopeX,1);
nSlopeY = size(slopeY,1);
nSlopes = nSlopeX + nSlopeY;
interactionMatrix = zeros( dm.nAct, nSlopeX+nSlopeY );

%% Make sure things are not saturated

% Record the current DM commands
curDMCommands = dm.cmdVector;
maxPixel = 0;
start(vid)
tempBackground=zeros(1024,1280,'uint8'); % Don't background subtract
for i=1:dm.nAct
    dm.cmdVector(i)=0.3;
    pause(0.1);
    [imageData,ts] = PWFSImageCaptureModulation(vid,tempBackground);
    curMaxPixel = max(max(imageData));
    if curMaxPixel > maxPixel
        maxPixel=curMaxPixel;
    end
    dm.cmdVector(1:dm.nAct)=curDMCommands;
end
disp(maxPixel)
stop(vid)

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
    [sPushX,sPushY] = PWFSGetSlopesModulation(vid,background);
    dm.cmdVector(i)=-pushPullValue;
    pause(pauseTime)
    [sPullX,sPullY] = PWFSGetSlopesModulation(vid,background);
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
condNum = 5.5; % Max Eigenvalue / Eigenvalue threshold
index = condVec > condNum;

% Plot the eigenvalues and the chosen conditioning parameter
f2 = figure('Name','Eigenvalues');
semilogy(eigenValues/eigenValues(1),'.')
line([1,97],[1./condNum,1./condNum],'Color','red')
xlabel('Eigenmodes')
ylabel('Normalized Eigenvalues')

% Threshold the eigenvalues
nThresholded = sum(index);
% nThresholded = 47;
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

%% Manual evaluation of the AO loop

% Declarations before evaluating the loop
dmIndex = 1:dm.nAct;
start(vid);
figure();

% Manual loop block
curVector = dm.cmdVector;
[sX,sY] = PWFSGetSlopesModulation(vid,background);
sVec = [sX;sY];
commandVec = curVector - loopGain*commandMatrix'*sVec;
dm.cmdVector(dmIndex) = commandVec;
[imageData,ts] = PWFSImageCaptureModulation(vid,background);
imageData = double(imageData);
imagesc(imageData)

% Include in block - Post correction slope standard deviation
% [sX,sY] = PWFSGetSlopesModulation(vid);
% fprintf('\nX std:%f Y std:%f\n',std(sX),std(sY))

% Reset the DM to flat and stop the camera
dm.Reset();
dm.cmdVector(1:dm.nAct) = curDMCommands;
stop(vid)

%% Create the loop
dmIndex = 1:dm.nAct;
start(vid);
pause(1);
figure();

nsteps = 100;
for i=1:nsteps

    curVector = dm.cmdVector;

    % Get slopes
    [sX,sY] = PWFSGetSlopesModulation(vid,background);
    sVec = [sX;sY];
    commandVec = curVector - loopGain*commandMatrix'*sVec;

    dm.cmdVector(dmIndex) = commandVec;

    pause(0.05)

    [imageData,ts] = PWFSImageCaptureModulation(vid,background);
    imageData = double(imageData);
    imagesc(imageData)

    disp(i)
end

dm.Reset();
dm.cmdVector(1:dm.nAct) = curDMCommands;
stop(vid)

%% How similar are the commands that would be applied?

% get a set of slopes
% sXTest = wfs.slopeX;
% sYTest = wfs.slopeY;
% sTest = [sXTest;sYTest];
%
% myCommands = commandMatrix'*sTest
% alpaoCommands = loop.commandMatrix'*sTest
%
% f9 = figure('Name','Compare DM commands');
% subplot(3,1,1)
% plot(1:dm.nAct, commandVec)
% xlabel('Actuator')
% ylabel('Poke magnitude')
%
% subplot(3,1,2)
% plot(1:dm.nAct, alpaoCommands)
% xlabel('Actuator')
% ylabel('Poke magnitude')
%
% subplot(3,1,3)
% for i=1:size(highVariationActuators,2)
%     hold on
%     thisHighVActuator = highVariationActuators(i);
%     plot([thisHighVActuator,thisHighVActuator],[-0.02,0.03],'Color',...
%          'Red','Linestyle','--','LineWidth',0.1)
% end
% hold on
% plot(1:dm.nAct, myCommands-alpaoCommands )
% hold off
% xlabel('Actuator')
% ylabel('Poke magnitude')

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
