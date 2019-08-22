% PWFSGetSlopes.m
% Get slopes

function [SxVec,SyVec] = PWFSGetSlopes(vid)

    % Get frame
    [imageData,ts] = PWFSImageCapture(vid);
    imageData = double(imageData);
    
    % Settings
    nPupil = 4;
    pupilExtractGeometry = 'circular'; % Geometry to extract pupils
    pupilRadius = 73; % Pixels
    pupilCol = [265,900,267,901];
    pupilRow = [200,197,833,828];
    pupilNames = ["Pupil 1","Pupil 2","Pupil 3","Pupil 4"];
    
    % Extract sub-images
    extractRadius = pupilRadius+5;
    I1 = imageData( pupilRow(1)-extractRadius:pupilRow(1)+extractRadius,...
                    pupilCol(1)-extractRadius:pupilCol(1)+extractRadius );
    I2 = imageData( pupilRow(2)-extractRadius:pupilRow(2)+extractRadius,...
                    pupilCol(2)-extractRadius:pupilCol(2)+extractRadius );
    I3 = imageData( pupilRow(3)-extractRadius:pupilRow(3)+extractRadius,...
                    pupilCol(3)-extractRadius:pupilCol(3)+extractRadius );
    I4 = imageData( pupilRow(4)-extractRadius:pupilRow(4)+extractRadius,...
                    pupilCol(4)-extractRadius:pupilCol(4)+extractRadius );
                
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

end