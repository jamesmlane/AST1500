%% measure_slopes.m
% Author - James Lane
% Capture an image of the PWFS pupils and measure the resulting slopes

%% Get the image
%
tiffObject = Tiff('./test_images/pyramid_test_2019-06-21-182902-0000.tif','r');
imageData = read(tiffObject);
imageData = sum(imageData,3);

%% Show the image

f1 = figure('Name','Pyramid WFS Image');
colormap gray
imagesc(imageData)

%% Define the location and size of the pupils

nPupil = 4;
pupilExtractGeometry = 'circular'; % Geometry to extract pupils
pupilRadius = 76; % Pixels
pupilCol = [330,973,333,978];
pupilRow = [195,190,835,825];
pupilNames = ["Pupil 1","Pupil 2","Pupil 3","Pupil 4"];

%% Draw circles around the pupil locations

f2 = figure('Name','Pyramid WFS Pupil Locations');
colormap gray
imagesc(imageData)
for i=1:nPupil
    hold on
    drawcircle(pupilCol(i),pupilRow(i),pupilRadius,pupilNames(i))
end
hold off
legend
    
%% Loop over each pupil and extract it

extractRadius = pupilRadius+5;

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
drawcircle(extractRadius,extractRadius,pupilRadius,pupilNames(1))
hold off
title(pupilNames(1))

subplot(2,2,2)
imagesc( I2 )
hold on
drawcircle(extractRadius,extractRadius,pupilRadius,pupilNames(2))
hold off
title(pupilNames(2))

subplot(2,2,3)
imagesc( I3 )
hold on
drawcircle(extractRadius,extractRadius,pupilRadius,pupilNames(3))
hold off
title(pupilNames(3))

subplot(2,2,4)
imagesc( I4 )
hold on
drawcircle(extractRadius,extractRadius,pupilRadius,pupilNames(4))
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

% First remove all of the negative values
I1(I1 < 0)=0;
I2(I2 < 0)=0;
I3(I3 < 0)=0;
I4(I4 < 0)=0;

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

function h = drawcircle(x,y,r,name)
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit, 'DisplayName',name);
end