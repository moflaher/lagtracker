function [lag,grid,time]=rungekutta(lag,grid,time)
%!==============================================================================|
%!  integrate particle position from x0 to xn using velocity fields at time     |
%!  t0 and time tn                                                              |
%!  npts:     number of particles					         |
%!  dti:   time step (usually dti)				                 |
%!  pdxn:     x position at time t0, returns x position at time tn              |
%!  pdyn:     y position at time t0, returns y position at time tn              |
%!  pdzn:     sigma position at time t0, returns sigma position at time tn      |
%!  pdznt:    z position at time tn                                             |
%!  host:     current set of elements containing particles                      |
%!  indomain: current status of particles, updated for particle position at tn  |
%!  sbound:   host on solide boundary                                           |
%!  u1/v1/w1: velocity field (u,v,omega) at time t0                             |
%!  u2/v2/w2: velocity field (u,v,omega) at time tn                             |
%!  hl:       bathymetry                                                        | 
%!  el1,el2:  free surface elevation at time t0 and tn                          | 
%!==============================================================================|
%!------------------------------------------------------------------------------|
%!  ns   : number of stages in explicit runga-kutta                             |
%!  chix : stage function evaluation for x-velocity                             | 
%!  pdx  : stage particle x-position                                            |
%!  ul   : stage u velocity                                                     |
%!  eps  : parameter defining depth of dry element                              |
%!  dmax : maximim sigma depth                                                  |
%!  a_rk : erk coefficients (a)                                                 |
%!  b_rk : erk coefficients (b)                                                 |
%!  c_rk : erk_coefficients (c)                                                 |
%!------------------------------------------------------------------------------|
 
	mstage=4;
	eps=1.0e-5;
	a_rk=[0,0.5,0.5,1];
	b_rk=[1/6,1/3,1/3,1/6];
	c_rk=[0,0.5,0.5,1]; 

%!--loop over rk stages 
lag.chix=0*lag.chix;
lag.chiy=0*lag.chiy;
lag.chiz=0*lag.chiz;
lag.chix(:,1) = lag.up(:);
lag.chiy(:,1) = lag.vp(:);
lag.chiz(:,1) = lag.wp(:)./(lag.hp+lag.ep);  


	for ns=2:mstage%mainfor
%!!particle position at stage n (x,y,sigma) 
		lag.xpt  = lag.xp(:)  + (a_rk(ns)*time.dti).*lag.chix(:,ns-1);
		lag.ypt  = lag.yp(:)  + (a_rk(ns)*time.dti).*lag.chiy(:,ns-1);
		lag.sigpt  = lag.sigp(:)  + (a_rk(ns)*time.dti).*lag.chiz(:,ns-1);

%!!adjust sigma position to stick to bottom
		lag.sigpt = max(lag.sigpt,-1);
%!!adjust sigma position to remain below free surface
		lag.sigpt = min(lag.sigpt,0);

%!!calculate velocity field for stage n using c_rk coefficients
		grid.uin  = ((1-c_rk(ns))*grid.u1 + c_rk(ns)*grid.u2)';
		grid.vin  = ((1-c_rk(ns))*grid.v1 + c_rk(ns)*grid.v2)'; 
		grid.win  = ((1-c_rk(ns))*grid.w1 + c_rk(ns)*grid.w2)'; 
		grid.ein  = (1-c_rk(ns))*grid.el1 + c_rk(ns)*grid.el2;  
		grid.uin(grid.nele+1,:)=0;    
		grid.vin(grid.nele+1,:)=0;
		grid.win(grid.nele+1,:)=0;

     
%!!evaluate velocity (u,v,w) at stage ns particle position
		[lag,grid]=newinterpolatev(lag,grid);

%!!evaluate el/h at stage ns particle position
		[lag,grid]=newinterpolateelh(lag,grid,0);
        %! evaluate the swim speed in sigma coordinates
		[lag]=calculate_swimspeed(lag);
 %!--above also add swimspeed to particles 
       lag.wp=lag.wp+lag.omega_swimspeed;
 
		lag.chix(:,ns) = lag.up(:);
		lag.chiy(:,ns) = lag.vp(:);
		lag.chiz(:,ns) = lag.wp(:)./(lag.hp+lag.ep);  

%!!limit vertical motion in very shallow water
		lag.chiz((lag.hp + lag.ep) < eps,ns)=0;
	end%mainfor

%!--sum stage contributions to get updated particle positions-------------------!
	lag.xpt  = lag.xp(:);
	lag.ypt  = lag.yp(:);
	lag.sigpt = lag.sigp(:);
	for ns=1:mstage
		lag.xpt = lag.xpt + time.dti*b_rk(ns)*lag.indomain(:).*lag.chix(:,ns);
		lag.ypt = lag.ypt + time.dti*b_rk(ns)*lag.indomain(:).*lag.chiy(:,ns);
		lag.sigpt = lag.sigpt + time.dti*b_rk(ns)*lag.indomain(:).*lag.chiz(:,ns);
	end


%!--evaluate temporary location
	
	[allfound,lag]=newfindquick(lag,grid);
	if(allfound==0)
		[lag, grid]=findfull(lag,grid,0);
	end

%determine if lines intersect turbines
    if lag.check_turbine_intersect
           lag=check_turbine_intersections(lag,grid);
    end

%!--update only particle still in water
			lag.xp  = lag.xp.*(1-lag.inwater)+lag.xpt.*lag.inwater;
			lag.yp  = lag.yp.*(1-lag.inwater)+lag.ypt.*lag.inwater;
			lag.sigp  = lag.sigp.*(1-lag.inwater)+lag.sigpt.*lag.inwater;
			lag.xpt  = lag.xp;
			lag.ypt  = lag.yp;
            

%!--adjust depth of updated particle positions----------------------------------!
	lag.sigp = max(lag.sigp(:),-1);	%stick to bottom
  	lag.sigp = min(lag.sigp(:),0);	%dont pierce free surface
	lag.sigpt  = lag.sigp;

%!--evaluate bathymetry and free surface height at updated particle position----!
	[lag,grid]=newinterpolatev(lag,grid);
	[lag,grid]=newinterpolateelh(lag,grid,0);
end


