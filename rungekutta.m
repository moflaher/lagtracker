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

if grid.diffusion
    lag.diffh=0*lag.diffh;
    lag.diffx=0*lag.diffx;
    lag.diffy=0*lag.diffy;
    lag.diffv=0*lag.diffv;
    lag.diffz=0*lag.diffz;

    lag.diffh(:,1) = lag.viscofhp(:);
    lag.diffx(:,1) = lag.viscofhx(:);
    lag.diffy(:,1) = lag.viscofhy(:);
    lag.diffv(:,1) = lag.khp(:)./(lag.hp+lag.ep); 
    lag.diffz(:,1) = lag.khz(:)./(lag.hp+lag.ep);  
end


	for ns=2:mstage%mainfor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particle position at stage n 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if grid.diffusion
            lag.xpt  = lag.xp(:)  + ((a_rk(ns)*time.dti).*(lag.chix(:,ns-1) + lag.diffx(:,ns-1))) + (sqrt(2*lag.diffh(:,ns-1))*lag.weiner((4*(time.iint-1))+ns));
		    lag.ypt  = lag.yp(:)  + ((a_rk(ns)*time.dti).*(lag.chiy(:,ns-1) + lag.diffy(:,ns-1))) + (sqrt(2*lag.diffh(:,ns-1))*lag.weiner((4*(time.iint-1))+ns));
		    lag.sigpt  = lag.sigp(:)  + ((a_rk(ns)*time.dti).*(lag.chiz(:,ns-1) + lag.diffz(:,ns-1))) + (sqrt(2*lag.diffv(:,ns-1))*lag.weiner((4*(time.iint-1))+ns));
        else       
            lag.xpt  = lag.xp(:)  + (a_rk(ns)*time.dti).*lag.chix(:,ns-1);
		    lag.ypt  = lag.yp(:)  + (a_rk(ns)*time.dti).*lag.chiy(:,ns-1);
		    lag.sigpt  = lag.sigp(:)  + (a_rk(ns)*time.dti).*lag.chiz(:,ns-1);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust sigma position to stick to bottom and remail below free surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		lag.sigpt = max(lag.sigpt,-1);
		lag.sigpt = min(lag.sigpt,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate velocity field for stage n using c_rk coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		grid.uin(1:grid.nele,:)  = ((1-c_rk(ns))*grid.u1 + c_rk(ns)*grid.u2);
		grid.vin(1:grid.nele,:)  = ((1-c_rk(ns))*grid.v1 + c_rk(ns)*grid.v2); 
		grid.win(1:grid.nele,:)  = ((1-c_rk(ns))*grid.w1 + c_rk(ns)*grid.w2); 
		grid.ein  = (1-c_rk(ns))*grid.el1 + c_rk(ns)*grid.el2;  

        if grid.diffusion
            grid.viscofhin  = ((1-c_rk(ns))*grid.viscofh1 + c_rk(ns)*grid.viscofh2);
		    grid.khin  = ((1-c_rk(ns))*grid.kh1 + c_rk(ns)*grid.kh2); 
        end	



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate velocity (u,v,w,el,h,dedtin) at stage ns particle position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[lag,grid]=newinterpolatev(lag,grid);
		[lag,grid]=newinterpolateelh(lag,grid,0);
        if grid.diffusion
            [lag,grid]=interpolate_diffusion(lag,grid,0);
        end	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate swimspeed, this is probably the wrong location for this which makes the units wrong...
% I think this should be change to be on line 39 with lag.chiz so that the velocity is only added once and in terms of m/s instead of d sig/ dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[lag]=calculate_swimspeed(lag);
        lag.wp=lag.wp+lag.omega_swimspeed;
 
		lag.chix(:,ns) = lag.up(:);
		lag.chiy(:,ns) = lag.vp(:);
		lag.chiz(:,ns) = lag.wp(:)./(lag.hp+lag.ep);  
        if grid.diffusion
            lag.diffh(:,ns) = lag.viscofhp(:);
            lag.diffx(:,ns) = lag.viscofhx(:);
            lag.diffy(:,ns) = lag.viscofhy(:);
            lag.diffv(:,ns) = lag.khp(:)./(lag.hp+lag.ep); 
            lag.diffz(:,ns) = lag.khz(:)./(lag.hp+lag.ep);  
        end

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
        if grid.diffusion
            lag.xpt = lag.xpt + (time.dti*b_rk(ns)*lag.indomain(:).*(lag.chix(:,ns) + lag.diffx(:,ns))) + (sqrt(2*lag.diffh(:,ns))*lag.weiner((4*(time.iint-1))+ns));
		    lag.ypt = lag.ypt + (time.dti*b_rk(ns)*lag.indomain(:).*(lag.chiy(:,ns) + lag.diffy(:,ns))) + (sqrt(2*lag.diffh(:,ns))*lag.weiner((4*(time.iint-1))+ns));
		    lag.sigpt = lag.sigpt + (time.dti*b_rk(ns)*lag.indomain(:).*(lag.chiz(:,ns) + lag.diffz(:,ns))) + (sqrt(2*lag.diffv(:,ns))*lag.weiner((4*(time.iint-1))+ns));
        else
		    lag.xpt = lag.xpt + time.dti*b_rk(ns)*lag.indomain(:).*lag.chix(:,ns);
		    lag.ypt = lag.ypt + time.dti*b_rk(ns)*lag.indomain(:).*lag.chiy(:,ns);
		    lag.sigpt = lag.sigpt + time.dti*b_rk(ns)*lag.indomain(:).*lag.chiz(:,ns);
        end
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
    if grid.diffusion
        [lag,grid]=interpolate_diffusion(lag,grid,0);
    end	

end


