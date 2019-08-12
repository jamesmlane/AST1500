function [images,timeStamp] = PWFSImageCapture_Modulation(vid)

% Flush the buffer
flushdata(vid);

% Wait until there's a frame available
while vid.FramesAvailable == 0
    pause(0.01)
end

% Read data
[images, timeStamp] = getdata(vid,1);

end
