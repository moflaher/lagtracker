function [lag,grid,time]=setlag(lag,grid,time)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:  lag,grid,time
% Return: lag,grid,time
%
%   Set the particle initial positions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lag.xp = lag.xpt(:);
	lag.yp = lag.ypt(:);
	lag.sigp  = lag.sigpt(:);    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine element containing each particle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	lag.ifound(:)=0;
	lag.indomain(:)=1;
	[lag grid]=findfull(lag, grid,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write initial particle positions and velocities to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lag.simstart=grid.time(1)+678942;
	time.itout=1;
	lag.x(:,time.itout)=lag.xpt(:)+ grid.vxmin;
	lag.y(:,time.itout)=lag.ypt(:)+ grid.vymin;
	lag.sig(:,time.itout)=lag.sigpt(:);
    [lag,grid]=newinterpolatev(lag,grid);
    [lag,grid]=newinterpolateelh(lag,grid,0);
    if grid.diffusion
        [lag,grid]=interpolate_diffusion(lag,grid,0);
    end	
	lag.u(:,time.itout)=lag.up;
	lag.v(:,time.itout)=lag.vp;
	lag.w(:,time.itout)=lag.wp;
	lag.z(:,time.itout)=lag.sigpt.*(lag.ep+lag.hp)+lag.ep;
	lag.time(time.itout)=time.starthour*time.instp;
    %lag.turbine_intersects=zeros(grid.nele,grid.nturbines);
    %lag.turbine_sigma=zeros(grid.nturbines,1);
    %lag.turbine_sigma=grid.zz(grid.turbine_sigmas(time.starthour+1,:));
end

	
	
           



