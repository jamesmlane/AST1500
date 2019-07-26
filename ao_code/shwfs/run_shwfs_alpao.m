% run_shwfs.m
% Author - James Lane
% Runs the Shack-Hartmann Wavefront Sensor by interfacing with the
% DM and HASO using ACE and implementing custom routines that will be
% eventually used for the pyramid.


% Initiating the system
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

loop.StartMonitoring();

% Initiate the loop
loop = aceloopLOOP();
loop.set('sWfs', wfs);
loop.set('sWfc', dm);
loop.Online();

% Build the IM
nAverage = 10;
nbPushPull= 2;
loop.BuildIM(0.1, nAverage, nbPushPull);
loop.BuildCM();
loop.BuildZ2C();
loop.SaveMatrix('C:\Users\Admin\Documents\MATLAB\MYMATRIX_20190529.mat');
loop.LoadMatrix('C:\Users\Admin\Documents\MATLAB\MYMATRIX_20180529.mat');

% Close the loop
loop.gain=0.3;
loop.Close();

loop.Open();
dm.Reset();

%%% IGNORE THIS PART
% Apply zeroes
dm.Reset();
dm.cmdVector=zeros(97,1);

% Save Flat command
loop.Open();
cv=dm.cmdVector;
save 'D:\Siqi\flat20181123.dat' cv -ascii
cv2= load('D:\Siqi\flat20181123.dat');
dm.cmdVector=cv2;

% Testing Zernike modes
oldZernikeVector = dm.zernikeVector;
test = zeros(1,30);
test(1)=0.8;
test(2)=1.3;
test(3) = 0.4;
test(4)=-0.01;
test(5)=-0.09;
test(6)=0.02;
test(7) = 0.05;
test(8) = 0.03;
test(9)=-0.06;
test(10)=-0.02;
test(11)=-0.1;
test(14)=0.05;
test(15)=0.06;
dm.zernikeVector = oldZernikeVector + test;
wfs.mode


%% Open loop & reset mirror
loop.Open();
bias = loop.sWfc.cmdVector; % Store current solution in bias variable
loop.sWfc.cmdVector = zeros(loop.sWfc.nAct, 1);

vid = videoinput('pointgrey', 1, 'F7_Mono12_1920x1200_Mode0');
src = getselectedsource(vid);

src.FrameRateMode = 'Manual';
src.GainMode = 'Manual';
src.ShutterMode = 'Manual';
src.ExposureMode = 'Off';
src.SharpnessMode = 'Off';

src.FrameRate = 14.9;
src.Shutter = 1.0;
src.Gain = 0;

snapshot = getsnapshot(vid);
imagesc(snapshot);

%% Set new bias and reference slopes
loop.sWfc.biasVector = bias;
% 4-point 5l/d
loop.SetReferenceZern([2.245 2.245 0 0 0 0 0 0 0 0 0]); % Close on 1um of focus
loop.SetReferenceZern([2.245 -2.245 0 0 0 0 0 0 0 0 0]); % Close on 1um of focus
loop.SetReferenceZern([-2.245 2.245 0 0 0 0 0 0 0 0 0]); % Close on 1um of focus
loop.SetReferenceZern([-2.245 -2.245 0 0 0 0 0 0 0 0 0]); % Close on 1um of focus
% 16-pt 5l/d
defocus = 0.0;
astig = 0.0;
lambda = 0.635;
loop.SetReferenceZern([4.904 0.975 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = getsnapshot(vid)*1.0;
loop.SetReferenceZern([4.157 2.778 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = im + getsnapshot(vid)*1.0;
loop.SetReferenceZern([2.778 4.157 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = im + getsnapshot(vid)*1.0;
loop.SetReferenceZern([0.975 4.904 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = im + getsnapshot(vid)*1.0;
loop.SetReferenceZern([-0.975 4.904 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = im + getsnapshot(vid)*1.0;
loop.SetReferenceZern([-2.778 4.157 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = im + getsnapshot(vid)*1.0;
loop.SetReferenceZern([-4.157 2.778 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = im + getsnapshot(vid)*1.0;
loop.SetReferenceZern([-4.904 0.975 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = im + getsnapshot(vid)*1.0;
loop.SetReferenceZern([-4.904 -0.975 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = im + getsnapshot(vid)*1.0;
loop.SetReferenceZern([-4.157 -2.778 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = im + getsnapshot(vid)*1.0;
loop.SetReferenceZern([-2.778 -4.157 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = im + getsnapshot(vid)*1.0;
loop.SetReferenceZern([-0.975 -4.904 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = im + getsnapshot(vid)*1.0;
loop.SetReferenceZern([0.975 -4.904 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = im + getsnapshot(vid)*1.0;
loop.SetReferenceZern([2.778 -4.157 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = im + getsnapshot(vid)*1.0;
loop.SetReferenceZern([4.157 -2.778 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = im + getsnapshot(vid)*1.0;
loop.SetReferenceZern([4.904 -0.975 defocus astig 0 0 0 0 0 0 0]*lambda);
pause(2);
im = im + getsnapshot(vid)*1.0;
imagesc(im);

imagesc(single(im)./single(noaber));

loop.Close();


%% Close the loop, and do slope calculation

loop.sWfc.biasVector = bias;
loop.Close();

% 16-pt 5l/d
lambda = 0.635;
im = zeros(1200,1920,'single');

defocus = 0.0;
astig = 1.0;
oastig = 0.0;

for i=1:16
    i=1;
    tip = cos(pi*(2*i-1)/16)*5;
    tilt = sin(pi*(2*i-1)/16)*5;
    loop.SetReferenceZern([tip tilt defocus astig oastig 0 0 0 0 0 0]*lambda);
    pause(2);
%    im = im + single(getsnapshot(vid))*1.0;
end
imagesc(im);

figure(7);
imagesc(single(im)./single(noaber));colorbar();hold on;
normed_im = single(im)./single(noaber);
normed_im = double(normed_im);
save 'D:\Siqi\Modulate\1l_astig_20180913\1l_astig_over_noaber.dat' normed_im -ascii


% Create a logical image of a circle with specified diameter, center, and image size.
% First create the image.
imageSizeX = 1920;
imageSizeY = 1200;
[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX = 729; centerY = 540;radius = 63;
circlePixels1 = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;% circlePixels is a 2D "logical" array.
pupil1 = normed_im.*circlePixels1;
center_pupil1 = pupil1(centerY-radius:centerY+radius, centerX-radius:centerX+radius);

centerX = 875; centerY = 540;radius = 63;
circlePixels2 = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;% circlePixels is a 2D "logical" array.
pupil2 = normed_im.*circlePixels2;
center_pupil2 = pupil2(centerY-radius:centerY+radius, centerX-radius:centerX+radius);

centerX = 729; centerY = 685;radius = 63;
circlePixels3 = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;% circlePixels is a 2D "logical" array.
pupil3 = normed_im.*circlePixels3;
center_pupil3 = pupil3(centerY-radius:centerY+radius, centerX-radius:centerX+radius);

centerX = 875; centerY = 685;radius = 63;
circlePixels4 = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2;% circlePixels is a 2D "logical" array.
pupil4 = normed_im.*circlePixels4;
center_pupil4 = pupil4(centerY-radius:centerY+radius, centerX-radius:centerX+radius);

sum_pupils = center_pupil1 + center_pupil2 + center_pupil3 + center_pupil4;
Sx = ((center_pupil1+center_pupil2)-(center_pupil3+center_pupil4))./sum_pupils;
Sx = double(Sx);
Sy = ((center_pupil1+center_pupil3)-(center_pupil2+center_pupil4))./sum_pupils;
Sy = double(Sy);

save 'D:\Siqi\Modulate\1l_astig_20180913\Sx_1l_astig.dat' Sx -ascii
save 'D:\Siqi\Modulate\1l_astig_20180913\Sy_1l_astig.dat' Sy -ascii

imagesc(Sy);
imagesc(center_pupil1);

Sx = load('D:\Siqi\Modulate\1l_astig_20180913\Sx_1l_astig.dat');
Sy = load('D:\Siqi\Modulate\1l_astig_20180913\Sy_1l_astig.dat');
normed_im = load('D:\Siqi\Modulate\1l_astig_20180913\1l_astig_over_noaber.dat');
%diff_x = diff(Sx, 1, 1);
histx=histogram(Sx,11);
hx_values = histx.Values;

%% Optical Gain
% Original slopes
defocus = 0.0;
astig = 1.0;
oastig = 0.0;
loop.SetReferenceZern([0 0 defocus astig oastig 0 0 0 0 0 0]*lambda);
Phase = dm.cmdMap;
imagesc(Phase);
diff_x = diff(Phase, 1, 1);
diff_y = diff(Phase, 1, 2);
figure(7);imagesc(diff_x);
figure(8);imagesc(diff_y);
%% Close the program

%
% IM=zeros(97,1);
% i=3;
% IM(i,1)=0.1;
% dm.cmdVector=cv2+IM;

% Shut-down the monitoring windows
loop.StopMonitoring();
loop.StartSaving();
loop.StopSaving();
loop.StartMonitoring();
loop.Open();
dm.Reset();

% Closing the windows
loop.StopMonitoring();
wfs.StopRtd();
wfs.StopSlopeRtd();
wfs.StopAlignmentRtd();
wfs.StopWavefrontRtd();

clear classes;
%clear;
