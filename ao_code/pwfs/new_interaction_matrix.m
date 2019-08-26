%% new_interaction_matrix.m
% Author - James Lane
% New code to generate an interaction matrix to try and get more modes

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