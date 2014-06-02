function [grid]=ncdread(grid,hour)


h=mod(hour,grid.size)+1;


grid.unc2=grid.unct(:,:,h);
grid.vnc2=grid.vnct(:,:,h);
grid.wnc2=grid.wnct(:,:,h);
grid.elnc2=grid.elnct(:,h);
