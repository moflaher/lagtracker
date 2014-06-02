function [lag,grid,time]=rungekutta(lag,grid,time)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: lag,grid,time
% Return: lag,grid,time
% 
%   Integrate all particle position using velocity fields grid.u1 and grid.u2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over rk stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lag.chix=0*lag.chix;
lag.chiy=0*lag.chiy;
lag.chiz=0*lag.chiz;
lag.chix(:,1) = lag.up(:);
lag.chiy(:,1) = lag.vp(:);
lag.chiz(:,1) = lag.wp(:)./(lag.hp+lag.ep);  



	for ns=2:mstage%mainfor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particle position at stage n 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		lag.xpt  = lag.xp(:)  + (a_rk(ns)*time.dti).*lag.chix(:,ns-1);
		lag.ypt  = lag.yp(:)  + (a_rk(ns)*time.dti).*lag.chiy(:,ns-1);
		lag.sigpt  = lag.sigp(:)  + (a_rk(ns)*time.dti).*lag.chiz(:,ns-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust sigma position to stick to bottom and remail below free surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		lag.sigpt = max(lag.sigpt,-1);
		lag.sigpt = min(lag.sigpt,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate velocity field for stage n using c_rk coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		grid.uin  = ((1-c_rk(ns))*grid.u1 + c_rk(ns)*grid.u2);
		grid.vin  = ((1-c_rk(ns))*grid.v1 + c_rk(ns)*grid.v2); 
		grid.win  = ((1-c_rk(ns))*grid.w1 + c_rk(ns)*grid.w2); 
		grid.ein  = (1-c_rk(ns))*grid.el1 + c_rk(ns)*grid.el2;  
		grid.uin(grid.nele+1,:)=0;    
		grid.vin(grid.nele+1,:)=0;
		grid.win(grid.nele+1,:)=0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate velocity (u,v,w,el,h,dedtin) at stage ns particle position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[lag,grid]=newinterpolatev(lag,grid);
		[lag,grid]=newinterpolateelh(lag,grid,0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate swimspeed, this is probably the wrong location for this which makes the units wrong...
% I think this should be change to be on line 39 with lag.chiz so that the velocity is only added once and in terms of m/s instead of d sig/ dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[lag]=calculate_swimspeed(lag);
        lag.wp=lag.wp+lag.omega_swimspeed;
 
		lag.chix(:,ns) = lag.up(:);
		lag.chiy(:,ns) = lag.vp(:);
		lag.chiz(:,ns) = lag.wp(:)./(lag.hp+lag.ep);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Limit vertical motion in very shallow water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		lag.chiz((lag.hp + lag.ep) < eps,ns)=0;
	end%mainfor



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sum stage contributions to get updated particle positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	lag.xpt  = lag.xp(:);
	lag.ypt  = lag.yp(:);
	lag.sigpt = lag.sigp(:);
	for ns=1:mstage
		lag.xpt = lag.xpt + time.dti*b_rk(ns)*lag.indomain(:).*lag.chix(:,ns);
		lag.ypt = lag.ypt + time.dti*b_rk(ns)*lag.indomain(:).*lag.chiy(:,ns);
		lag.sigpt = lag.sigpt + time.dti*b_rk(ns)*lag.indomain(:).*lag.chiz(:,ns);
	end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate temporary location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	[allfound,lag]=newfindquick(lag,grid);
	if(allfound==0)
		[lag, grid]=findfull(lag,grid,0);
	end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine if lines intersect turbines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
    if lag.check_turbine_intersect
           lag=check_turbine_intersections(lag,grid);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update only the particles still in water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			lag.xp  = lag.xp.*(1-lag.inwater)+lag.xpt.*lag.inwater;
			lag.yp  = lag.yp.*(1-lag.inwater)+lag.ypt.*lag.inwater;
			lag.sigp  = lag.sigp.*(1-lag.inwater)+lag.sigpt.*lag.inwater;
			lag.xpt  = lag.xp;
			lag.ypt  = lag.yp;
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust depth of updated particle positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	lag.sigp = max(lag.sigp(:),-1);	%stick to bottom
  	lag.sigp = min(lag.sigp(:),0);	%dont pierce free surface
	lag.sigpt  = lag.sigp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate velocity (u,v,w,el,h,dedtin) at updated particle position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	[lag,grid]=newinterpolatev(lag,grid);
	[lag,grid]=newinterpolateelh(lag,grid,0);

end


