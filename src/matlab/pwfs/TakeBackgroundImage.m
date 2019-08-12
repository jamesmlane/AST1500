%% TakeBackgroundImage.m
% Take control of the laser and take a background image

function [background] = TakeBackgroundImage(vid)

  fig=figure();
  lasersource=actxcontrol('APTLASER.APTLaserCtrl.1',[20 20 600 400],fig); %connect to ActivX controller
  lasersource.StartCtrl %start control

  SN = 86874274; % Set the Serial Number
  set(lasersource,'HWSerialNum', SN);

  lasersource.Identify; %identify the device
  pause(2); % waiting for the GUI to load up;

  lasersource.DisableOutput; %shutting off the laser

  pause(0.01);


  %taking the background image
  flushdata(vid);
  pause(0.04)
  wait(vid,0.05,'logging')
  [background, timeStamp] = getdata(vid,1);
  pause(0.1);

  %Setting up the laser & turning it back on
  lasersource.SetGUIDispUnits(2); %setting units to mW (3 is dBm)
  lasersource.SetPowerSetpoint(0.5); %0 to 4.9mW
  lasersource.EnableOutput(); %turn on the laser

  lasersource.StopCtrl; %stop controlling the laser
  delete(lasersource);
  clear lasersource;

end
