% make_cm.m
% Author - James Lane
% Make the command matrix for the Shack Hartmann Wavefront Sensor

% Define how many pushes and pulls there will be, and what their magnitudes
% will be.
pauseTime = 0.1;      % Time to pause between command and slope acquisition
nAverage = 5;         % Number of frames to average for each push-pull
nbPushPull = 5;       % Number of push-pull
pushPullValue = 0.12; % Magnitude of push-pull

% Setup the interaction matrix array
nSlopeX = size(wfs.slopeX,1);
nSlopeY = size(wfs.slopeY,1);
nSlopes = nSlopeX + nSlopeY;
interactionMatrix = zeros( dm.nAct, nSlopeX+nSlopeY );

% Loop over the number of actuators
for i=1:dm.nAct
  
  % Setup the array that will hold the push-pull samples
  sDiff = zeros(nSlopes,1)
  
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
    
    % Calculate the difference and append
    sDiffX = sPushX-sPullX;
    sDiffY = sPushY-sPullY;
    sDiff = sDiff + [sDiffX,sDiffY];
    
  end
  
  % Take the average of the slope vectors
  sVec = sDiff / nbPushPull / pushPullValue

  % Append this column into the interaction matrix
  interactionMatrix(i,:) = sVec;
end
