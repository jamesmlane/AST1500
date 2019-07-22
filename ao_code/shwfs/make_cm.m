%% make_cm.m
% Author - James Lane
% Make the command matrix for the Shack Hartmann Wavefront Sensor

%% Initiating the system
% getenv('ACEROOT')
setenv( 'ACEROOT', 'C:\Users\Admin\Documents\Alpao\AlpaoCoreEngine')
addpath( fullfile(getenv('ACEROOT'), 'matlab') );
acecsStartup();
userStartup
dm.Reset();

% Switch on the monitoring windows
wfs.StartRtd();
wfs.StartSlopeRtd();
wfs.StartAlignmentRtd();
wfs.StartWavefrontRtd();
dm.StartMonitoring();
wfs.sCam.dit = 0.06; % Set camera exposure time in ms

% Initiate the loop
loop.StartMonitoring();
loop = aceloopLOOP();
loop.set('sWfs', wfs);
loop.set('sWfc', dm);
loop.Online();

%% Definitions
% Number of pushes and pulls, and their magnitudes
pauseTime = 0.1;      % Time to pause between command and slope acquisition
nAverage = 5;         % Number of frames to average for each push-pull
nbPushPull = 5;       % Number of push-pull
pushPullValue = 0.12; % Magnitude of push-pull

% Setup the interaction matrix array
nSlopeX = size(wfs.slopeX,1);
nSlopeY = size(wfs.slopeY,1);
nSlopes = nSlopeX + nSlopeY;
interactionMatrix = zeros( dm.nAct, nSlopeX+nSlopeY );

%% Make the ALPAO command matrix
loop.BuildIM(0.1, nAverage, nbPushPull);
loop.BuildCM();
loop.BuildZ2C();

%% Create our interaction matrix
wbfig = waitbar(0,'Creating interaction matrix...');

% Loop over the number of actuators
for i=1:dm.nAct
  
  % Setup the array that will hold the push-pull samples
  sDiff = zeros(nSlopes,1);
  
  % For each actuator loop over the number of push-pull pairs
  for j=1:nbPushPull
    
    % Acquire slopes
    dm.cmdVector(i)=pushPullValue;
    pause(pauseTime)
    sPushX = wfs.slopeX;
    sPushY = wfs.slopeY;
    dm.cmdVector(i)=-pushPullValue;
    pause(pauseTime)
    sPullX = wfs.slopeX;
    sPullY = wfs.slopeY;
    dm.cmdVector(i)=0;
    
    % Calculate the difference and append
    sDiffX = sPushX-sPullX;
    sDiffY = sPushY-sPullY;
    sDiff = sDiff + [sDiffX;sDiffY];
    
  end
  
  % Take the average of the slope vectors
  sVec = sDiff/(2*nbPushPull)/pushPullValue;

  % Append this column into the interaction matrix
  interactionMatrix(i,:) = sVec;
  
  waitbar(i/dm.nAct)
end
close(wbfig)

%% Plots of the interaction matrices

f1 = figure('Name','My Interaction Matrix');
imagesc(interactionMatrix);
colorbar

f2 = figure('Name','ALPAO Interaction Matrix');
imagesc(loop.interactionMatrix);
colorbar

%% Perform SVD and pseudoinverse
[U,S,V] = svd(interactionMatrix,'econ');
eigenValues = diag(S);
iS = diag(1./eigenValues);

subplot(1,2,2)
semilogy(eigenValues,'.')
xlabel('Eigen modes')
ylabel('Eigen values')

% The last 4 eigen values are filtered out



%% Shut it down
dm.Reset();
dm.StopMonitoring();
wfs.StopRtd();
wfs.StopSlopeRtd();
wfs.StopAlignmentRtd();
wfs.StopWavefrontRtd();
loop.StopMonitoring();