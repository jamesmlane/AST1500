%% CameraTriggerLoop.m
% Flush data from the camera buffer

function CameraTriggerLoop(vid)

flushdata(vid,'trigger');

end
