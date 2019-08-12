%% measure_slopes.m
% Author - James Lane
% Capture an image of the PWFS pupils and measure the resulting slopes

%% Set local paths
addpath( genpath('../../src/matlab/pwfs') )

%% Start the video object
vid = videoinput('pointgrey', 1);
get(vid)
shutter = 0.01;
flushdata(vid);% clears all frames from buffer
src = getselectedsource(vid);
vid.FramesPerTrigger = 1;
vid.TriggerRepeat = inf;
triggerconfig(vid,'manual');
start(vid);

%% Get the image
%
% tiffObject = Tiff('./test_images/pyramid_test_2019-06-21-182902-0000.tif','r');
% imageData = read(tiffObject);
% imageData = sum(imageDat)

flushdata(vid);
[imageData,ts] = PWFSImageCapture(vid);
imageData = double(imageData);

%% Show the image

f1 = figure('Name','Pyramid WFS Image');
colormap gray
imagesc(imageData)

%% Define the location and size of the pupils

nPupil = 4;
pupilExtractGeometry = 'circular'; % Geometry to extract pupils
pupilRadius = 73; % Pixels
pupilCol = [265,900,267,901];
pupilRow = [200,197,833,828];
pupilNames = ["Pupil 1","Pupil 2","Pupil 3","Pupil 4"];

%% Draw circles around the pupil locations

f2 = figure('Name','Pyramid WFS Pupil Locations');
colormap gray
imagesc(imageData)
for i=1:nPupil
    hold on
    draw_circle(pupilCol(i),pupilRow(i),pupilRadius,pupilNames(i))
end
hold off
legend

%% Loop over each pupil and extract it

extractRadius = pupilRadius+10;

I1 = imageData( pupilRow(1)-extractRadius:pupilRow(1)+extractRadius,...
                pupilCol(1)-extractRadius:pupilCol(1)+extractRadius );
I2 = imageData( pupilRow(2)-extractRadius:pupilRow(2)+extractRadius,...
                pupilCol(2)-extractRadius:pupilCol(2)+extractRadius );
I3 = imageData( pupilRow(3)-extractRadius:pupilRow(3)+extractRadius,...
                pupilCol(3)-extractRadius:pupilCol(3)+extractRadius );
I4 = imageData( pupilRow(4)-extractRadius:pupilRow(4)+extractRadius,...
                pupilCol(4)-extractRadius:pupilCol(4)+extractRadius );

f3 = figure('Name','Pyramid WFS Pupil Extraction')
colormap gray

subplot(2,2,1)
imagesc( I1 )
hold on
draw_circle(extractRadius,extractRadius,pupilRadius,pupilNames(1))
hold off
title(pupilNames(1))

subplot(2,2,2)
imagesc( I2 )
hold on
draw_circle(extractRadius,extractRadius,pupilRadius,pupilNames(2))
hold off
title(pupilNames(2))

subplot(2,2,3)
imagesc( I3 )
hold on
draw_circle(extractRadius,extractRadius,pupilRadius,pupilNames(3))
hold off
title(pupilNames(3))

subplot(2,2,4)
imagesc( I4 )
hold on
draw_circle(extractRadius,extractRadius,pupilRadius,pupilNames(4))
hold off
title(pupilNames(4))

%% Slope Calculations

% Define the valid pixel map
[xExtract,yExtract] = meshgrid(1:1+2*extractRadius,1:1+2*extractRadius);
xExtract = xExtract - extractRadius;
yExtract = yExtract - extractRadius;
rExtract = sqrt( xExtract.^2 + yExtract.^2 );
validPixelMap = rExtract < pupilRadius;
nValidPixels = size(find(validPixelMap));

% Now set all pixels outside of the pupil mask to 0
I1(~validPixelMap)=0;
I2(~validPixelMap)=0;
I3(~validPixelMap)=0;
I4(~validPixelMap)=0;

% INorm = I1+I2+I3+I4; % Normalize by total flux
INorm = 0.25*mean( I1(validPixelMap)+I2(validPixelMap)+I3(validPixelMap)+...
              I4(validPixelMap) )*ones(size(validPixelMap));

% I think I've made sure to define this properly. See Figure 2.
% It's at least internally consistent, but I don't think it's consistent
% in an absolute sense with Verinaud or OOMAO geometry
SyMap = ( (I1+I2) - (I3+I4) );
SxMap = ( (I1+I3) - (I2+I4) );

SyVec = SyMap(validPixelMap);
SxVec = SxMap(validPixelMap);

f4 = figure('Name','Slope Maps');

subplot(1,2,1)
colormap parula
imagesc( SxMap )
colorbar
title('X Slopes')

subplot(1,2,2)
colormap parula
imagesc( SyMap )
colorbar
title('Y Slopes')

%% Functions

function h = draw_circle(x,y,r,name)
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit, 'DisplayName',name);
end
