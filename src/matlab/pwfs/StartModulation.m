function [E727,Controller] = StartModulation(frequency, amplitude, xscale, phaseshift)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BeginWithThisSampleWhenUsing_E727.m
%
% This sample demonstrates how to
% - load the PI MATLAB Driver GCS2
% - connect to the E727
% - initialize the E727
% - set the stage connected to E727
% - move the E727
% - configure and use the data recorder
% - disconnect the E727
% - unload the PI MATLAB Driver GCS2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load PI MATLAB Driver GCS2

if (strfind(evalc('ver'), 'Windows XP'))
    if (~exist('C:\Documents and Settings\All Users\PI\PI_MATLAB_Driver_GCS2','dir'))
        error('The PI MATLAB Driver GCS2 was not found on your system. Probably it is not installed. Please run PI_MATLAB_Driver_GCS2_Setup.exe to install the driver.');
    else
        addpath('C:\Documents and Settings\All Users\PI\PI_MATLAB_Driver_GCS2');
    end
elseif (strfind(evalc('ver'), 'Windows'))
    if (~exist('C:\Users\Public\PI\PI_MATLAB_Driver_GCS2','dir'))
        error('The PI MATLAB Driver GCS2 was not found on your system. Probably it is not installed. Please run PI_MATLAB_Driver_GCS2_Setup.exe to install the driver.');
    else
        addpath('C:\Users\Public\PI\PI_MATLAB_Driver_GCS2');
    end
end


if(~exist('Controller','var'))
    Controller = PI_GCS_Controller();
end
if(~isa(Controller,'PI_GCS_Controller'))
    Controller = PI_GCS_Controller();
end


%% List connected USB and TCP/IP controller

% devicesUsb = Controller.EnumerateUSBAsArray();
% disp([num2str(length(devicesUsb)), ' PI controllers are connected to your PC via USB']);
% disp(devicesUsb);
%
% devicesTcpIp = Controller.EnumerateTCPIPDevicesAsArray();
% disp([num2str(length(devicesTcpIp)), ' PI controllers were found on the network']);
% disp(devicesTcpIp);


%% Parameters
% You MUST EDIT AND ACITVATE the parameters to make your system run properly.
% 1. Set the correct stage type. AN INCORRECT STAGE TYPE  CAN DAMAGE YOUR STAGE!
% 2. Activate the connection type you want to use.
% 3. Set the connection settings.


% Set the correct stage type. A WRONG STAGE TYPE CAN DAMAGE YOUR STAGE!
stageType = 'S-330.4SH';

% Connection settings
use_RS232_Connection    = false;
use_USB_Connection      = true;
use_TCPIP_Connection    = false;


if (use_RS232_Connection)
    comPort = 1;          % Look at the Device Manager to get the rigth COM Port.
    baudRate = 115200;    % Look at the manual to get the rigth bau rate for your controller.
end

if (use_USB_Connection)
    controllerSerialNumber = '0118001659'; %'Enter valid serial number here, e.g. "123456789"';    % Use "devicesUsb = Controller.EnumerateUSB('')" to get all PI controller connected to you PC.
                                                                                    % Or look at the label of the case of your controller
end

if (use_TCPIP_Connection)
    ip = '192.168.2.109';  % Use "devicesTcpIp = Controller.EnumerateTCPIPDevices('')" to get all PI controller available on the network.
    port = 50000;           % Is 50000 for almost all PI controllers
end


%% Start connection

boolE727connected = false;

if (exist('E727','var'))
    if (E727.IsConnected)
        boolE727connected = true;
    end
end


if (~boolE727connected)
    if (use_RS232_Connection)
        E727 = Controller.ConnectRS232(comPort, baudRate);
    end

    if (use_USB_Connection)
        E727 = Controller.ConnectUSB(controllerSerialNumber);
    end

    if (use_TCPIP_Connection)
        E727 = Controller.ConnectTCPIP(ip, port);
    end
end


% query controller identification
E727.qIDN()

% initialize controller
E727 = E727.InitializeController();


%% Configuration and referencing

% query controller axes
availableaxes = E727.qSAI_ALL; % changed from E727.qSAI_ALLasArray;
if(isempty(availableaxes))
	error('No axes available');
end
%axisname = availableaxes{1};
axis_one = availableaxes{1};
axis_two = availableaxes{2};
axis_thr = availableaxes{3};


% Auto Zero all available axes
E727.ATZ(availableaxes)

% wait until all axes are set to zero: all([1 1 1]) >> 1
while(0 == all(E727.qATZ()))
    pause(0.1);
end

% Switch on servo
E727.SVO(axis_one, 1); % 1 for servo on
E727.SVO(axis_two, 1);


%% Basic controller functions
dMin_one = E727.qTMN(axis_one); % 0 to 5000
dMax_one = E727.qTMX(axis_one);
dMin_two = E727.qTMN(axis_two); % 0 to 5000
dMax_two = E727.qTMX(axis_two);
dMin_thr = E727.qTMN(axis_thr); % 0 to 100
dMax_thr = E727.qTMX(axis_thr);

%position_one = rand(1)*(dMax-dMin)+dMin;
position_one =0.5*(dMax_one-dMin_one)+dMin_one;
% amplitude of sine waves are 5 with phase difference of lambda/4;
% offset by 2.5 units to correct starting position:
position_two = 0.5*(dMax_two-dMin_two)+dMin_two+2.5;
E727.GcsCommandset('POS?');
pause(1);
disp(E727.GcsGetAnswer())
% SVA to do open-loop
E727.MOV(axis_one, position_one);
E727.GcsCommandset('POS?');
pause(1);
disp(E727.GcsGetAnswer())
E727.MOV(axis_two, position_two);

E727.qPOS(axis_one)
E727.qPOS(axis_two)

%{
%%Hysteresis step code (Remeber to set SVO to 0 and change MOV to SVA%%

posinitial=[position_one,E727.qPOS(axis_one),position_two,E727.qPOS(axis_two)];

dlmwrite('Hysteresis.txt',posinitial,'-append','delimiter',' ','roffset',0);

pause(1);


step=linspace(100,5000,50);

for i=1:50

    E727.SVA(axis_one,position_one+step(i));
%E727.qPOS(axis_one)
    pause(1);
    E727.SVA(axis_two,position_two+step(i));
%E727.qPOS(axis_two)
    pause(1);

    motion1=[position_one+step(i),E727.qPOS(axis_one),position_two+step(i),E727.qPOS(axis_two)];

    dlmwrite('Hysteresis.txt',motion1,'-append','delimiter',' ','roffset',0);
    pause(1);
end

for i=1:49

    E727.SVA(axis_one,position_one +step(end-i));
%E727.qPOS(axis_one)
    pause(1);
    E727.SVA(axis_two,position_two+step(end-i));
%E727.qPOS(axis_two)
    pause(1);

    motion2=[position_one+step(end-i),E727.qPOS(axis_one),position_two+step(end-i),E727.qPOS(axis_two)];

    dlmwrite('Hysteresis.txt',motion2,'-append','delimiter',' ','roffset',0);
    pause(1);
end


E727.SVA(axis_one,position_one);
E727.SVA(axis_two,position_two);

posfinal=[position_one,E727.qPOS(axis_one),position_two,E727.qPOS(axis_two)];

dlmwrite('Hysteresis.txt',posfinal,'-append','delimiter',' ','roffset',0);

stop
%}

disp(position_one);
disp(dMin_one);
disp(dMax_one);

disp(position_two);
disp(dMin_two);
disp(dMax_two);

% wait for motion to stop
while(0 ~= E727.IsMoving(axis_one) || 0 ~= E727.IsMoving(axis_two))
    pause(0.001); % seconds
end

E727.qPOS(axis_one)
E727.qPOS(axis_two)
E727.qPOS(axis_thr)
%% Send an receive GCS-Commands directly
% This is especially usefull if an specific GCS-Command isn't implemented
% in the MATLAB-Driver yet. GCS-Commands are described in the user manual
% of the controller.

% The following 4 lines do (nearly) the same as "E727.qPOS()" does:
E727.GcsCommandset('POS?');
pause(1);
disp(E727.GcsGetAnswer())
errorNumber = E727.GetError()
if (0 ~= errorNumber), error(E727.TranslateError(errorNumber)); end;


%% Call Function
relativeMove = 1;   %in �m
dataSource = 2;     %Current Position of axis,
% PI_DataRecorderStep(E727, axisname, relativeMove, dataSource);

%E727.SVO(axis_one,1);
%pause(1);
%E727.SVO(axis_two,1);
pause(1);

E727.DRC(1,axis_one,1)% target position
E727.DRC(2,axis_one,2)% current position
E727.DRC(3,axis_two,1)
E727.DRC(4,axis_two,2)

figure(1);


%E727.STE(axisname,1)

%while(E727.IsMoving(axisname))
%    pause(0.1);
%end


%Setting up data recorder

E727.CCL(1,'advanced'); %setting command permissions to level 1

E727.DRC(transpose([1 2 5 7]),axis_one,transpose([1 2 3 27])); %definining which data sets to record (x axis)
E727.DRC(transpose([3 4 6 8]),axis_two,transpose([1 2 3 27])); %definining which data sets to record (y axis)


% define sine wave for wave table 1 (WAV 1):
    % segment length = 1000
    % amplitude = 5
    % offset= 2500 (mid-stroke)
    % wavelength = 1000
    % start point = 0
    % curve center point = 500

%Phasemap=importdata('C:\Users\Admin\Desktop\LUKE\BodeData_1lamdaD_betterPID.txt');




SUT=5e-5; %seconds

%Modulation parameters
freq=frequency;%370;%Hz
amp=amplitude;%3000; %diameter of modulation (urad)
%circular corrections
ygravcorrect=1.0;%0.9;
xcorrect=xscale;%cos(45*pi/180)
phasex=phaseshift;%30;
phasey=0;%-20;


wavelength=round(1/(freq*SUT)); %1kHz=20
offset=2500-(amp/2); %lowest value possible, not midpoint
seglength=wavelength;
curvecenterpoint=wavelength/2;
startpointx=0;
startpointy=wavelength/4; %offset one axis


%Generating custom wave tables [WAV_PNT(tableid,tablestartpoint(must be 1),#ofpts,append/destroyparam(0=destroy,1=append),arrayofwavepoints)]
E727.WAV_PNT(1,1,wavelength,0,SinWaveGen(freq,amp,2500,xcorrect,phasex));
E727.WAV_PNT(2,1,wavelength,0,SinWaveGen(freq,amp,2500,1,90));



%Generating the wave tables
%WaveX=['WAV 1 X SIN_P ',int2str(seglength),' ',int2str(xcorrect*amplitude),' ',int2str(offset + 0.5*(1-xcorrect)*amplitude),' ',int2str(wavelength),' ',int2str(startpointx-(phasex/(360*frequency*SUT))), ' ',int2str(curvecenterpoint)];
%WaveY=['WAV 2 X SIN_P ',int2str(seglength),' ',int2str(ygravcorrect*amplitude),' ',int2str(offset + 0.5*(1-ygravcorrect)*amplitude),' ',int2str(wavelength),' ',int2str(startpointy-(phasey/(360*frequency*SUT))), ' ',int2str(curvecenterpoint)];


%E727.GcsCommandset('WTR 0 1 1') %increase the wave 1 output duration by factor of 3
%E727.GcsCommandset('WTR 0 0.5 2') %increase the wave 2 output duration by factor of 3





%E727.GcsCommandset(WaveX)
%disp(E727.GetError())
%pause(1);
%E727.GcsCommandset(WaveY)
%disp(E727.GetError())

pulsehigh=zeros(4);
pulsehigh(end)=seglength;
for i = 1:4
pulsehigh(end-i)=seglength-i;
end

%Setting up digital output trigger
E727.TWC() %clear all output trigger settings so all points are set to "LOW"

E727.TWS(1,seglength,1);
%E727.TWS(transpose([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]),transpose([seglength,seglength-1,seglength-2,seglength-3,seglength-4,seglength-5,seglength-6,seglength-7,seglength-8,seglength-9,seglength-10,seglength-11,seglength-12,seglength-13,seglength-14,seglength-15,seglength-16,seglength-17,seglength-18,seglength-19]),transpose([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])); %set digital output to be high when wave 1 reaches the last point in the wave table and low everywhere else

E727.CTO(1,3,4);%connect to output trigger 3(9 for short pulse, 4 for sUT pulse)

E727.CTO(1,7,0);
%moving the stage to the start point of the wave table
%Xwavedata=E727.qGWD(1,seglength,[1]); %get Xwave data
%Ywavedata=E727.qGWD(1,seglength,[2]); % get Ywave data

%E727.MOV(axis_one, Xwavedata(1));%move axes to 1st data point
%E727.MOV(axis_two, Ywavedata(1));

%pause(0.2);

% connect wave generator (axis 1 & 2) to wave table 1 & 2
E727.GcsCommandset('WSL 1 1')
disp(E727.GetError())
E727.GcsCommandset('WSL 2 2')
disp(E727.GetError())






% start output of wave generators 1 & 2
pause(2)

E727.GcsCommandset('WGO 1 1')
E727.GcsCommandset('WGO 2 1')
disp(E727.GetError())

end
