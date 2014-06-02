function [grid]=ncdread(grid,nh)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: grid,nh
% Return: grid
% 
%  Updates velocity and elevation fields to time nh.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncid=netcdf.open(grid.ncfile,'NC_NOWRITE');

grid.unc2 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[0 0 nh],[grid.nele grid.siglay 1]);
grid.vnc2 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[0 0 nh],[grid.nele grid.siglay 1]);
grid.wnc2 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ww'),[0 0 nh],[grid.nele grid.siglay 1]);
grid.elnc2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zeta'),[0 nh],[grid.node 1]);

netcdf.close(ncid);
