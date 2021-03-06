%% run_shwfs.m
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
wfs.sCam.dit = 0.06; % Set camera exposure time in ms

% Initiate the loop
loop.StartMonitoring();
loop = aceloopLOOP();
loop.set('sWfs', wfs);
loop.set('sWfc', dm);
loop.Online();

% Monitor the DM
dm.StopMonitoring();
dm.StartMonitoring();

%% Definitions
% Number of pushes and pulls, and their magnitudes, loop gain
pauseTime = 0.1;      % Time to pause between command and slope acquisition
nAverage = 5;         % Number of frames to average for each push-pull
nbPushPull = 5;       % Number of push-pull
pushPullValue = 0.12; % Magnitude of push-pull
loopGain = 0.3;       % Loop gain

% Setup the interaction matrix array
nSlopeX = size(wfs.slopeX,1);
nSlopeY = size(wfs.slopeY,1);
nSlopes = nSlopeX + nSlopeY;
interactionMatrix = zeros( dm.nAct, nSlopeX+nSlopeY );

%% Make the ALPAO command matrix
% loop.BuildIM(0.1, nAverage, nbPushPull);
% loop.BuildCM();
% loop.BuildZ2C();

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
xlabel('Lenslet Slope Measurement')
ylabel('Actuator')
imagesc(interactionMatrix);
colorbar

% f2 = figure('Name','ALPAO Interaction Matrix');
% xlabel('Lenslet Slope Measurement')
% ylabel('Actuator')
% imagesc(loop.interactionMatrix);
% colorbar

%% Perform SVD and pseudoinverse
[U,S,V] = svd(interactionMatrix);
eigenValues = diag(S);
iS = diag(1./eigenValues);

% Decide on the number of eigenvalues to threshold
condVec = eigenValues(1)./eigenValues;
condNum = 4; % Max Eigenvalue / Eigenvalue threshold
index = condVec > condNum;
nThresholded = sum(index);
% nThresholded = 24
fprintf('%i / 97 modes removed\n',nThresholded)

% Plot the eigenvalues
f3 = figure('Name','Eigenvalues')
semilogy(eigenValues/eigenValues(1),'.')
line([1,97],[1./condNum,1./condNum],'Color','red')
xlabel('Eigenmodes')
ylabel('Normalized Eigenvalues')

% Remake the inverse eigenvalue matrix with thresholding
iSth = diag(1./eigenValues(1:end-nThresholded));
[nS,nC] = size(interactionMatrix);
iSth(nC,nS) = 0;

% Compute the command matrix
commandMatrix = V*iSth*U';

%% Plot the command matrices

f5 = figure('Name','My Command Matrix');
imagesc(commandMatrix);
xlabel('Actuator')
ylabel('Lenslet Slope Measurement')
colorbar
caxis([-70,70])

% f6 = figure('Name','ALPAO Command Matrix');
% imagesc(loop.commandMatrix);
% xlabel('Actuator')
% ylabel('Lenslet Slope Measurement')
% colorbar
% caxis([-70,70])

% f7 = figure('Name','Command Matrix Difference');
% diffCommandMatrix = loop.commandMatrix - commandMatrix;
% imagesc(diffCommandMatrix);
% xlabel('Actuator')
% ylabel('Lenslet Slope Measurement')
% colorbar
% caxis([-10,10])

% f8 = figure('Name','Standard Deviation of Actuator Response');
% stdDevCommandMatrixColumns = zeros(dm.nAct,1);
% for i=1:dm.nAct
%     diffCommandMatrixColumn = diffCommandMatrix(:,i);
%     stdDevCommandMatrixColumns(i) = std(diffCommandMatrixColumn);
% end
% plot(1:dm.nAct,stdDevCommandMatrixColumns,'.')
% xlabel('Actuator')
% ylabel('Slope Response Standard Deviation')

% hVAIndex = 1:dm.nAct;
% highVariationActuators = hVAIndex(stdDevCommandMatrixColumns > 2);

%% Create the loop

% First do a single application of the vector

% Get slopes
sX = wfs.slopeX;
sY = wfs.slopeY;
sVec = [sX;sY];

% Make DM commands
commandVec = commandMatrix'*sVec;
dmIndex = 1:dm.nAct;

% Apply the DM commands
dm.cmdVector(dmIndex) = -commandVec;


nsteps = 5000;
for i=1:nsteps
   
    curVector = dm.cmdVector;
    
    % Get slopes
    sX = wfs.slopeX;
    sY = wfs.slopeY;
    sVec = [sX;sY];
    
    commandVec = curVector - loopGain*commandMatrix'*sVec;
    
    dm.cmdVector(dmIndex) = commandVec;
    
    pause(0.01)
    
    disp(i)
end



%% How similar are the commands that would be applied?

% get a set of slopes
% sXTest = wfs.slopeX;
% sYTest = wfs.slopeY;
% sTest = [sXTest;sYTest];
%
% myCommands = commandMatrix'*sTest
% alpaoCommands = loop.commandMatrix'*sTest
%
f9 = figure('Name','Compare DM commands');
subplot(3,1,1)
plot(1:dm.nAct, commandVec)
xlabel('Actuator')
ylabel('Poke magnitude')
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
dm.Reset();
dm.StopMonitoring();
wfs.StopRtd();
wfs.StopSlopeRtd();
wfs.StopAlignmentRtd();
wfs.StopWavefrontRtd();
loop.StopMonitoring();
