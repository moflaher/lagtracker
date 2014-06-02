function [lag,grid,time]=setlag(lag,grid,time)
  
%!----shift coordinates to model coordinate system------------------------------|
%	lag.xp = lag.xpt(:) - grid.vxmin;
%	lag.yp = lag.ypt(:) - grid.vymin;
	lag.xp = lag.xpt(:);
	lag.yp = lag.ypt(:);
    
%!----determine element containing each particle--------------------------------|
	lag.ifound(:)=0;
	lag.indomain(:)=1;
	[lag grid]=findfull(lag, grid,1);
%   done in full find
%   lag.indomain(lag.ifound==0)=0;

	lag.sigp  = lag.sigpt(:);
	
 
%!----store current position array in previous time level array-----------------|  
 
	lag.xplst = lag.xp;
	lag.yplst = lag.yp;
	lag.sigplst = lag.sigp;

%!----count number of particles outside domain----------------------------------|
	lag.np_out = lag.npts - sum(lag.indomain);

%!------------------------------------------------------------------------------|
%!  print statistics on lagrangian tracking to output                           |
%!------------------------------------------------------------------------------|
            
% 	if(lag.np_out ~= 0)
%      		['# pts outside domain  :' int2str(lag.np_out) ]    
%   	else
% 		['# pts outside domain  :  none']    
%  	end

%!------------------------------------------------------------------------------|
%!  open up output files                                                        |
%!------------------------------------------------------------------------------|

%!--write initial particle positions and velocities to file

	time.itout=1;
	lag.x(:,time.itout)=lag.xpt(:)+ grid.vxmin;
	lag.y(:,time.itout)=lag.ypt(:)+ grid.vymin;
	lag.sig(:,time.itout)=lag.sigpt(:);
    [lag,grid]=newinterpolatev(lag,grid);
    [lag,grid]=newinterpolateelh(lag,grid,0);
	lag.u(:,time.itout)=lag.up;
	lag.v(:,time.itout)=lag.vp;
	lag.w(:,time.itout)=lag.wp;
	lag.z(:,time.itout)=lag.sigpt.*(lag.ep+lag.hp)+lag.ep;
	lag.time(time.itout)=time.starthour*time.instp;
    lag.turbine_intersects=zeros(grid.nele,grid.nturbines);
    lag.turbine_sigma=zeros(grid.nturbines,1);
    lag.turbine_sigma=grid.zz(grid.turbine_sigmas(time.starthour+1,:));
end

	
	
           



