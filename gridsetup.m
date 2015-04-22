function [grid]=gridsetup(time,set)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: time,set
% Return: grid
%
%   Loads data from ncfile and initializes arrays that will be used.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncid=netcdf.open(set.fvcompath,'NC_NOWRITE');
grid=gridsetupdef;

grid.ncfile=set.fvcompath;
if set.diffusion
    grid.diffusion=true;
else
    grid.diffusion=false;
end
if set.twoDflag==true
    grid.istwoD=true;
end

grid.vx=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'x'));
grid.vy=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'y'));
grid.nv=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'nv')));
grid.h=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'));
grid.hele=(grid.h(grid.nv(:,1))+grid.h(grid.nv(:,2))+grid.h(grid.nv(:,1)))/3;
temp=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'siglev'));
[nodes lsig]=size(temp);
grid.zz=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'siglay'),[0 0],[1 lsig-1])';
grid.z=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'siglev'),[0 0],[1 lsig])';
grid.nbe=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'nbe')));
grid.a1u=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'a1u'));
grid.a2u=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'a2u'));
grid.awx=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'awx'));
grid.awy=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'awy'));
grid.aw0=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'aw0'));
grid.nele=length(grid.nv);
grid.node=length(grid.h);
grid.siglay=length(grid.zz);
grid.siglev=length(grid.z);
grid.isonb=zeros(grid.nele,1);	
grid.isbce=zeros(grid.nele,1);
grid.xc=zeros(grid.nele,1);
grid.yc=zeros(grid.nele,1);
grid.uin=zeros(grid.nele+1,grid.siglay);
grid.vin=zeros(grid.nele+1,grid.siglay);		
grid.win=zeros(grid.nele+1,grid.siglay);
grid.hin=zeros(grid.node,1);
grid.hin=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'));
grid.ein=zeros(grid.node,1);
grid.dedtin=zeros(grid.node,1);

grid.unc1 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[0 0 set.start],[grid.nele grid.siglay 1]);
grid.vnc1 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[0 0 set.start],[grid.nele grid.siglay 1]);
if grid.istwoD
    grid.wnc1 = zeros(grid.nele,grid.siglay);
else
    grid.wnc1 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ww'),[0 0 set.start],[grid.nele grid.siglay 1]);
end
grid.elnc1=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zeta'),[0 set.start],[grid.node 1]);
grid.time=netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'));

grid.uin=grid.unc1;    
grid.vin=grid.vnc1;    
grid.win=grid.wnc1;  
grid.ein=grid.elnc1; 
grid.uin(grid.nele+1,:)=0;    
grid.vin(grid.nele+1,:)=0;
grid.win(grid.nele+1,:)=0;  

grid.unc2=zeros(grid.nele,grid.siglay);
grid.vnc2=zeros(grid.nele,grid.siglay);
grid.wnc2=zeros(grid.nele,grid.siglay);
grid.elnc2=zeros(grid.node,1);

grid.u1=zeros(grid.nele,grid.siglay);
grid.v1=zeros(grid.nele,grid.siglay);
grid.w1=zeros(grid.nele,grid.siglay);
grid.el1=zeros(grid.node,1);
grid.u2=zeros(grid.nele,grid.siglay);
grid.v2=zeros(grid.nele,grid.siglay);
grid.w2=zeros(grid.nele,grid.siglay);
grid.el2=zeros(grid.node,1);



if grid.diffusion
    grid.viscofhnc1 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'viscofh'),[0 0 set.start],[grid.node grid.siglay 1]);
    grid.khnc1 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'kh'),[0 0 set.start],[grid.node grid.siglev 1]);
    grid.viscofhnc2=zeros(grid.node,grid.siglay);
    grid.khnc2=zeros(grid.node,grid.siglev);
    
    grid.viscofhin=grid.viscofhnc1;    
    grid.khin=grid.khnc1;         
    
    grid.viscofh1=zeros(grid.node,grid.siglay);
    grid.kh1=zeros(grid.node,grid.siglev);
    grid.viscofh2=zeros(grid.node,grid.siglay);
    grid.kh2=zeros(grid.node,grid.siglev);
end


grid.vxmax=0;
grid.vxmin=0;
grid.vymax=0;
grid.vymin=0; 

%--compute differnce and intra-sigma levels
    grid.dzz(1:(grid.siglay-1))=grid.zz(1:(grid.siglay-1))-grid.zz(2:grid.siglay);
    grid.dz(1:(grid.siglev-1))=grid.z(1:(grid.siglev-1))-grid.z(2:grid.siglev);

    turb=[];%ncfvcom{'turbdrag'}(:);
    grid.turbines=[];%find(turb(1,grid.siglay,:)>0);
    grid.nturbines=length(grid.turbines);
    if ~isempty(grid.turbines)
        a=turb(:,:,grid.turbines)>0;
        grid.turbine_sigmas=squeeze(11-sum(a,2));
    end


netcdf.close(ncid);
end
  



