
function [images,timeStamp] = PWFSImageCapture(vid)

flushdata(vid);
trigger(vid);
    
while islogging(vid) == 1
    % disp('Capturing images...')
    pause(0.01)
end
    
[images, timeStamp] = getdata(vid,vid.FramesAvailable);


 %connecting



%vid=videoinput('pointgrey',1);

%triggerconfig(vid2, 'hardware', 'risingEdge', 'externalTriggerMode1-Source2');%setting trigger for 2nd camera
%configuring camera settings
%vid.ROIPosition = [830 585 96 96];



%src.FrameRateMode = 'Manual';
%src.GainMode = 'Manual'; 
%src.ShutterMode = 'Manual';
%src.ExposureMode = 'Off';
%src.SharpnessMode = 'Off';



%src.FrameRate = framerate; %87=minimum but this might vary based on the shutter
%src.Shutter = shutter; %0.01=minimum
%src.Gain = 0;
%src.Exposure=0.0083;

%vid.LoggingMode = 'disk'; %saving images to disk instead of memory
%logfile = VideoWriter('C:\Users\Admin\Desktop\LUKE\ModulationFlatImage.mj2', 'Motion JPEG 2000'); %setting the file type
%vid.DiskLogger = logfile; %creating the image logging object


%triggerconfig(vid,'manual') %specifying trigger type






%vid.TriggerFrameDelay = 0; %number of initially skipped frames

%vid.TriggerRepeat = Nimages-1;

%src.TriggerDelayMode='Off';

%triggerconfig(vid,'hardware','fallingEdge','externalTriggerMode0-Source0');

%images=ones(524,640,1,19);



%for i=1:19
%pause(1);

%src.Shutter=shutter(i);
%pause(0.5);



% if isrunning(vid) == 1 %confirming that the video feed is runnning
%     disp('Video feed enabled.')
% else 
%     disp('No video feed active.')
%     
% end

 %begin image acquisition
%{
while vid.TriggersExecuted == 0
    disp('waiting for trigger...')
    wait(vid,5);
end
%}

% pausing for image capture


% if vid.FramesAcquired == 0 %checking that trigger was successful 
%     disp('Trigger failure: 0 frames acquired.')
% else 
%     disp([int2str(vid.FramesAcquired),'','frames acquired.'])
% end

 
 %stop video feed

%images(:,:,:,i)=getdata(vid,vid.FramesAvailable);

%end

% [images, timeStamp] = getdata(vid,vid.FramesAvailable); %Move frames from buffer into MATLAB

%waittime = src.FrameRate * (vid.FramesPerTrigger + vid.TriggerFrameDelay) + 0; 
%wait(vid, waittime); %including wait time for buffering

%savedframes=vid.FramesAvailable;
%vid.DiskLoggerFrameCount;


% disp([int2str( vid.DiskLoggerFrameCount),' ','frames saved to memory.']) %checking number of saved images 




%figure; %display image
%imagesc(frames(:,:,:,1)); %frames saved as (height x width x colorband x frame)



%testing the centroid function
%[Xvals,Yvals]=Centroid2d(frames(:,:,:,1),350,7);

%Centroidmap=ones(1200,1920,'uint16');
  
%Centroidmap(Xvals,Yvals)=0;


%figure;
%imagesc((Centroidmap.*frames(:,:,:,1)));
 




end



