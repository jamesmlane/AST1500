function [cleanframes,timeStamp] = PWFSImageCapture_Modulation(vid,background)

flushdata(vid);

%pause(0.01)

wait(vid,0.05,'logging')

% while vid.FramesAvailable == 0
%     pause(0.01)
% end

% wait_time = 0;
% while islogging(vid) == 1 | wait_time > 5*exposure_time
%     % disp('Capturing images...')
%     pause(0.01)
%     wait_time = wait_time + 0.01;
% end
    
[images, timeStamp] = getdata(vid,1);

cleanframes=images(:,:,1,1)-background;


end