%% DrawCircle.m
% Utility function to draw a circle on a matlab plot

function h = draw_circle(x,y,r,name)
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit, 'DisplayName',name);
end
