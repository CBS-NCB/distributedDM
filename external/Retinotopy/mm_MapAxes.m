function mm_MapAxes(gca, nx, ny, mm_step)
%MA adopted from spatial_axis.m on 20150302
%   put mm on spatial axis of tensor images, according to mmperpix
%   nx and ny are the dimensions of the map
%   mm_step is the step size in millimeter
if nargin<4
    mm_step  = .5; % millimeter
end

%%
mmperpix = 0.0063; % 1/157.8 millimeter, see screen info or VDAQ
space = mmperpix*(1:max(nx,ny));
xx = 1;
while xx*mm_step < mmperpix*ny
    x_space_tick(xx) = find(space > xx*mm_step, 1);
    xx= xx+1;
end

yy=1;
while yy*mm_step < mmperpix*nx
    y_space_tick(yy) = find(space > yy*mm_step,1);
    yy= yy+1;
end

set(gca,'xtick', x_space_tick, 'xticklabel', round(10*space(x_space_tick))/10);
set(gca,'ytick', y_space_tick, 'yticklabel', round(10*space(y_space_tick))/10);
xlabel('mm')
ylabel('mm')

end