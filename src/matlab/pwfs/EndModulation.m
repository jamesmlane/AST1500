function EndModulation(E727,Controller)

% stop output of wave generators 1 & 2
E727.GcsCommandset('WGO 1 0')
E727.GcsCommandset('WGO 2 0')
disp(E727.GetError())

pause(0.5)

%% If you want to close the connection
E727.CloseConnection;

%% if you want to unload the dll and destroy the class object
Controller.Destroy;
clear Controller;
clear E727;
end
