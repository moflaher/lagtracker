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
	eps=.01;
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
lag.chiz(:,1) = lag.wp(:);  

if grid.diffusion
    ff=lag.savediffusionfudgefactor;
    lag.diffh=0*lag.diffh;
    lag.diffx=0*lag.diffx;
    lag.diffy=0*lag.diffy;
    lag.diffv=0*lag.diffv;
    lag.diffz=0*lag.diffz;

    lag.diffh(:,1) = lag.viscofhp(:)/ff;
    lag.diffx(:,1) = lag.viscofhx(:)/ff;
    lag.diffy(:,1) = lag.viscofhy(:)/ff;
    lag.diffv(:,1) = lag.khp(:)/ff; 
    lag.diffz(:,1) = lag.khz(:)/ff;  



end


	for ns=2:mstage%mainfor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particle position at stage n 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if grid.diffusion
            lag.xpt  = lag.xp(:)  + ((a_rk(ns)*time.dti).*(lag.chix(:,ns-1) + lag.diffx(:,ns-1))) + (sqrt(2*lag.diffh(:,ns-1))*lag.wiener((4*(time.iint-1))+ns-1));
		    lag.ypt  = lag.yp(:)  + ((a_rk(ns)*time.dti).*(lag.chiy(:,ns-1) + lag.diffy(:,ns-1))) + (sqrt(2*lag.diffh(:,ns-1))*lag.wiener((4*(time.iint-1))+ns-1));
		    lag.zpt  = lag.zp(:)  + ((a_rk(ns)*time.dti).*(lag.chiz(:,ns-1) + lag.diffz(:,ns-1))) + (sqrt(2*lag.diffv(:,ns-1))*lag.wiener((4*(time.iint-1))+ns-1));   
        else       
            lag.xpt  = lag.xp(:)  + (a_rk(ns)*time.dti).*lag.chix(:,ns-1);
		    lag.ypt  = lag.yp(:)  + (a_rk(ns)*time.dti).*lag.chiy(:,ns-1);
		    lag.zpt  = lag.zp(:)  + (a_rk(ns)*time.dti).*lag.chiz(:,ns-1);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust z position to stick to bottom and remail below free surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[lag,grid]=newinterpolateelh(lag,grid,1);
		lag.zpt = min(lag.zpt,lag.ep);
        
        %if particle is within 1 cm of bottom stop movement
        cutoff=.01;
        lag.inwater=((lag.zpt+lag.hp)>cutoff);

		lag.zpt = max(lag.zpt,-lag.hp);

        lag.sigpt=lag.zpt(:)./(-1*(lag.hp+lag.ep));

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
		[lag,grid]=newinterpolatev(lag,grid,0);
		[lag,grid]=newinterpolateelh(lag,grid,0);
        if grid.diffusion
            [lag,grid]=interpolate_diffusion(lag,grid,0);
        end	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate swimspeed (let the calculate_swinspeed function in incase it becomes more complex later)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[lag]=calculate_swimspeed(lag);
        lag.wp=lag.wp+lag.w_swimspeed;
 
		lag.chix(:,ns) = lag.up(:);
		lag.chiy(:,ns) = lag.vp(:);
		lag.chiz(:,ns) = lag.wp(:);  
        if grid.diffusion
            lag.diffh(:,ns) = lag.viscofhp(:)/ff;
            lag.diffx(:,ns) = lag.viscofhx(:)/ff;
            lag.diffy(:,ns) = lag.viscofhy(:)/ff;
            lag.diffv(:,ns) = lag.khp(:)/ff; 
            lag.diffz(:,ns) = lag.khz(:)/ff;  
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
	lag.zpt  = lag.zp(:);
	for ns=1:mstage
        if grid.diffusion
            lag.xpt = lag.xpt + (time.dti*b_rk(ns)*lag.indomain(:).*(lag.chix(:,ns) + lag.diffx(:,ns))) + lag.indomain(:).*(sqrt(2*lag.diffh(:,ns))*lag.wiener((4*(time.iint-1))+ns));
		    lag.ypt = lag.ypt + (time.dti*b_rk(ns)*lag.indomain(:).*(lag.chiy(:,ns) + lag.diffy(:,ns))) + lag.indomain(:).*(sqrt(2*lag.diffh(:,ns))*lag.wiener((4*(time.iint-1))+ns));
		    lag.zpt = lag.zpt + (time.dti*b_rk(ns)*lag.indomain(:).*(lag.chiz(:,ns) + lag.diffz(:,ns))) + lag.indomain(:).*(sqrt(2*lag.diffv(:,ns))*lag.wiener((4*(time.iint-1))+ns));
        else
		    lag.xpt = lag.xpt + time.dti*b_rk(ns)*lag.indomain(:).*lag.chix(:,ns);
		    lag.ypt = lag.ypt + time.dti*b_rk(ns)*lag.indomain(:).*lag.chiy(:,ns);
		    lag.zpt = lag.zpt + time.dti*b_rk(ns)*lag.indomain(:).*lag.chiz(:,ns);
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
			lag.zp  = lag.zp.*(1-lag.inwater)+lag.zpt.*lag.inwater;
			lag.xpt  = lag.xp;
			lag.ypt  = lag.yp;
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust depth of updated particle positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
    lag.zp = min(lag.zp,lag.ep); %dont pierce free surface
    lag.zp = max(lag.zp,-lag.hp); %stick to bottom
	lag.zpt  = lag.zp;

    lag.sigpt=lag.zpt(:)./(-1*(lag.hp+lag.ep));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate velocity (u,v,w,el,h,dedtin) at updated particle position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	[lag,grid]=newinterpolatev(lag,grid,1);
	[lag,grid]=newinterpolateelh(lag,grid,0);
    if grid.diffusion
        [lag,grid]=interpolate_diffusion(lag,grid,0);
    end	

end


