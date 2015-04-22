function [grid]=ncdread(grid,nh)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: grid,nh
% Return: grid
% 
%  Updates velocity and elevation fields to time nh.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncid=netcdf.open(grid.ncfile,'NC_NOWRITE');

if grid.diffusion

    grid.unc2 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[0 0 nh],[grid.nele grid.siglay 1]);
    grid.vnc2 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[0 0 nh],[grid.nele grid.siglay 1]);
    if grid.istwoD
        grid.wnc2 = zeros(grid.nele,grid.siglay);   
    else
        grid.wnc2 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ww'),[0 0 nh],[grid.nele grid.siglay 1]);
    end  
    grid.elnc2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zeta'),[0 nh],[grid.node 1]);
    grid.viscofhnc2 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'viscofh'),[0 0 nh],[grid.node grid.siglay 1]);
    grid.khnc2 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'kh'),[0 0 nh],[grid.node grid.siglev 1]);

else

    grid.unc2 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[0 0 nh],[grid.nele grid.siglay 1]);
    grid.vnc2 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[0 0 nh],[grid.nele grid.siglay 1]);
    if grid.istwoD
        grid.wnc2 = zeros(grid.nele,grid.siglay);   
    else
        grid.wnc2 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ww'),[0 0 nh],[grid.nele grid.siglay 1]);
    end    
    grid.elnc2=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zeta'),[0 nh],[grid.node 1]);

end

netcdf.close(ncid);
