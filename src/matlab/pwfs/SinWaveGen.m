%% SinWaveGen.m
% Generates a custom sine wave to run the 

function output = SinWaveGen(frequency,diameter,offset,ampscale,phaseshift)

  Npoints=1/(frequency*5e-5);

  angles=linspace(0,360-(frequency*360*5e-5),Npoints)*(pi/180);

  phase=phaseshift*(pi/180);

  output=(0.5*diameter*ampscale*(-1)*sin(angles+phase))+offset;

end
