function [grid]=gridsetup(time,set)

%nele               number of elements 
%node               number of nodes

%--------------------------grid metrics---------------------------------------------
%vxmin,vymin,vxmax,vymax	
%xc(:)               x-coord at face center 
%yc(:)               y-coord at face center
%vx(:)               x-coord at grid point
%vy(:)               y-coord at grid point

%----------------node, boundary condition, and control volume-----------------------
%nv(:,:)             node numbering for elements
%nbe(:,:)            indices of element neighbors
%ntve(:)         
%isonb(:)            node marker = 0,1,2   
%isbce(:)     
%nbve(:,:)
%nbvt(:,:)

%----------------1-d arrays for the sigma coordinate -------------------------------
%z(:)                    sigma coordinate value 
%zz(:)                   intra level sigma value
%dz(:)                   delta-sigma value
%dzz(:)                  delta of intra level sigma 


%---------------2-d flow variable arrays at nodes----------------------------------
%h(:)            bathymetric depth   
%el2(:)           current surface elevation
%el1(:)           surface elevation at previous time step

%---------------- internal mode   arrays-(element based)----------------------------
%u2(:,:)         x-velocity
%v2(:,:)         y-velocity
%w2(:,:)         omega-velocity
%u1(:,:)        x-velocity from previous timestep
%v1(:,:)        y-velocity from previous timestep
%w1(:,:)        omega velocity from previous time step

%------------shape coefficient arrays and control volume metrics--------------------
%a1u(:,:)      
%a2u(:,:)     
%awx(:,:)   
%awy(:,:)  
%aw0(:,:)



ncfvcom=netcdf(set.fvcompath,'nowrite');

grid=gridsetupdef;

grid.size=set.numberoftimesteps;

grid.vx=ncfvcom{'x'}(:);
grid.vy=ncfvcom{'y'}(:);
grid.nv=ncfvcom{'nv'}(:)';
grid.h=ncfvcom{'h'}(:);
grid.z=ncfvcom{'siglev'}(:);
grid.zz=ncfvcom{'siglay'}(:);
grid.nbe=ncfvcom{'nbe'}(:)';
grid.a1u=ncfvcom{'a1u'}(:)';
grid.a2u=ncfvcom{'a2u'}(:)';
grid.awx=ncfvcom{'awx'}(:)';
grid.awy=ncfvcom{'awy'}(:)';
grid.aw0=ncfvcom{'aw0'}(:)';

grid.nele=length(grid.nv);
grid.node=length(grid.h);
grid.siglay=length(ncfvcom{'siglay'}(:));%kbm1
grid.siglev=length(ncfvcom{'siglev'}(:));%kb

grid.isonb=zeros(grid.nele,1);	
grid.isbce=zeros(grid.nele,1);
grid.xc=zeros(grid.nele,1);
grid.yc=zeros(grid.nele,1);
grid.uin=zeros(grid.nele+1,grid.siglay);
grid.vin=zeros(grid.nele+1,grid.siglay);		
grid.win=zeros(grid.nele+1,grid.siglev);

grid.hin=zeros(grid.node,1);
grid.hin=ncfvcom{'h'}(:);
grid.ein=zeros(grid.node,1);
grid.dedtin=zeros(grid.node,1);
%grid.ntve=zeros(grid.node,1);

grid.unct=zeros(grid.siglay,grid.nele,grid.size);
grid.vnct=zeros(grid.siglay,grid.nele,grid.size);
grid.wnct=zeros(grid.siglev,grid.nele,grid.size);

%reshape matrices
unct=ncfvcom{'u'}(1:grid.size,1:grid.siglay,1:grid.nele);
vnct=ncfvcom{'v'}(1:grid.size,1:grid.siglay,1:grid.nele);
wnct=ncfvcom{'w'}(1:grid.size,1:grid.siglev,1:grid.nele);
for sigi=1:grid.siglay
    grid.unct(sigi,1:grid.nele,1:grid.size)=squeeze(unct(:,sigi,:))';
    grid.vnct(sigi,1:grid.nele,1:grid.size)=squeeze(vnct(:,sigi,:))';
end
for sigi=1:grid.siglev
    grid.wnct(sigi,1:grid.nele,1:grid.size)=squeeze(wnct(:,sigi,:))';
end

grid.elnct=(ncfvcom{'zeta'}(1:grid.size,1:grid.node))';

grid.unc1=grid.unct(:,:,time.starthour+1);
grid.vnc1=grid.vnct(:,:,time.starthour+1);
grid.wnc1=grid.wnct(:,:,time.starthour+1);
grid.elnc1=grid.elnct(:,time.starthour+1);

grid.uin(1:grid.nele,:)=grid.unc1';    
grid.vin(1:grid.nele,:)=grid.vnc1';    
grid.win(1:grid.nele,:)=grid.wnc1';    
grid.ein=grid.elnc1;    

grid.unc2=zeros(grid.siglay,grid.nele);
grid.vnc2=zeros(grid.siglay,grid.nele);
grid.wnc2=zeros(grid.siglev,grid.nele);
grid.elnc2=zeros(grid.node,1);

grid.u1=zeros(grid.siglay,grid.nele);
grid.v1=zeros(grid.siglay,grid.nele);
grid.w1=zeros(grid.siglev,grid.nele);
grid.el1=zeros(grid.node,1);
grid.u2=zeros(grid.siglay,grid.nele);
grid.v2=zeros(grid.siglay,grid.nele);
grid.w2=zeros(grid.siglev,grid.nele);
grid.el2=zeros(grid.node,1);

grid.vxmax=0;
grid.vxmin=0;
grid.vymax=0;
grid.vymin=0; 

%--compute differnce and intra-sigma levels
    grid.dzz(1:(grid.siglay-1))=grid.zz(1:(grid.siglay-1))-grid.zz(2:grid.siglay);
    grid.dz(1:(grid.siglev-1))=grid.z(1:(grid.siglev-1))-grid.z(2:grid.siglev);

    turb=ncfvcom{'turbdrag'}(:);
    grid.turbines=find(turb(1,grid.siglay,:)>0);
    grid.nturbines=length(grid.turbines);
    if ~isempty(grid.turbines)
        a=turb(:,:,grid.turbines)>0;
        grid.turbine_sigmas=squeeze(11-sum(a,2));
    end

close(ncfvcom);
end
  



