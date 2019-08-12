%% test_frame_grab.m
% Stress test the image capture routines to make sure the system doesn't crash
% This routine is somewhat obsolete now that we are using modulation

% Setup the video object
vid = videoinput('pointgrey', 1);

get(vid)
flushdata(vid); % clears all frames from buffer

shutter = 0.01;
src = getselectedsource(vid);
vid.FramesPerTrigger = 1;
vid.TriggerRepeat = inf;
triggerconfig(vid,'manual');

start(vid);

% Loop many times over the
for i=1:10000
    trigger(vid);

    while islogging(vid) == 1
        disp('Capturing images...')
        pause(1.5*shutter)
    end

    [images, timeStamp] = getdata(vid,vid.FramesAvailable);

    disp(i)
end

stop(vid);

delete(vid); %closing the connection
clear vid;
