function [lag,grid,time]=lagupdate(lag,grid,time)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: lag,grid,time
% Return: lag,grid,time
%
%   Main body of lagtracker.
%   This function loops over all times and calls rungekutta to calculate particle postions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update particle positions, calculate scalar fields and particle velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	grid.u1=grid.unc1;
	grid.v1=grid.vnc1;
	grid.w1=grid.wnc1;
    grid.el1=grid.elnc1;
    if grid.diffusion
        grid.viscofh1=grid.viscofhnc1;
        grid.kh1=grid.khnc1;
    end
    if grid.wave
        grid.u1Wave=grid.unc1Wave;
        grid.v1Wave=grid.vnc1Wave;
    end
    if grid.wind
        grid.u1Wind=grid.unc1Wind;
        grid.v1Wind=grid.vnc1Wind;
    end
	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over the tracking period  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for nh=time.starthour+1:time.finishhour%mainfor
		

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read velocity and elevation fields for next timestep from fvcompath
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[grid]=ncdread(grid,nh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop to move particle. Use linearly interpolated based on set.internalstep.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		for it=1:time.i2%main2

			['On step: ' int2str(time.iint+1) '/' int2str(time.i2*time.trackingtime) ]

			time.iint=time.iint+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time interpolation of physical fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        		tmp2=it/time.i2;
        		tmp1=(it-time.i2)/-time.i2;

        		grid.u2=tmp1*grid.unc1+tmp2*grid.unc2;
                grid.v2=tmp1*grid.vnc1+tmp2*grid.vnc2;
        		grid.w2=tmp1*grid.wnc1+tmp2*grid.wnc2;
        		grid.el2=tmp1*grid.elnc1+tmp2*grid.elnc2;
                grid.dedtin=(grid.el2-grid.el1)/time.dti;
                if grid.diffusion
                    grid.viscofh2=tmp1*grid.viscofhnc1+tmp2*grid.viscofhnc2;
                    grid.kh2=tmp1*grid.khnc1+tmp2*grid.khnc2;
                end
                if grid.wave
        		    tmp2=it/time.i2;
        		    tmp1=(it-time.i2)/-time.i2;
                    grid.u2Wave=tmp1*grid.unc1Wave+tmp2*grid.unc2Wave;
                    grid.v2Wave=tmp1*grid.vnc1Wave+tmp2*grid.vnc2Wave;
                end
                if grid.wind
        		    tmp2=it/time.i2;
        		    tmp1=(it-time.i2)/-time.i2;
                    grid.u2Wind=tmp1*grid.unc1Wind+tmp2*grid.unc2Wind;
                    grid.v2Wind=tmp1*grid.vnc1Wind+tmp2*grid.vnc2Wind;
                end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Move particle.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			[lag, grid, time]=rungekutta(lag,grid,time);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save data to lag at set.outputstep.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        		if(mod(time.iint,time.int2) == 0)%if3  
				    time.itout=time.itout+1;
                    array=lag.indomain==1;
                    lag.x(array,time.itout)=lag.xp(array) + grid.vxmin;
               		lag.y(array,time.itout)=lag.yp(array) + grid.vymin;
				    lag.sig(array,time.itout)=lag.sigp(array);
				    lag.z(array,time.itout)=lag.zp(array);
				    lag.u(array,time.itout)=lag.up(array);
               		lag.v(array,time.itout)=lag.vp(array);
				    lag.w(array,time.itout)=lag.wp(array);
                    lag.h(array,time.itout)=lag.hp(array);
                    array=lag.indomain==0;
                    lag.x(array,time.itout)=NaN;
               		lag.y(array,time.itout)=NaN;
				    lag.sig(array,time.itout)=NaN;
				    lag.z(array,time.itout)=NaN;
				    lag.u(array,time.itout)=0;
               		lag.v(array,time.itout)=0;
				    lag.w(array,time.itout)=0;
                    lag.h(array,time.itout)=NaN;
				    lag.time(time.itout)=(time.iint*time.dti)+(time.starthour*time.instp);
        		end%if3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time step update of velocity fields within move particle loop.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       		grid.u1=grid.u2;
			grid.v1=grid.v2;
			grid.w1=grid.w2;
			grid.el1=grid.el2;
            if grid.diffusion
                grid.viscofh1=grid.viscofh2;
                grid.kh1=grid.kh2;
            end
            if grid.wave
                grid.u1Wave=grid.u2Wave;
                grid.v1Wave=grid.v2Wave;
            end
            if grid.wind
                grid.u1Wind=grid.u2Wind;
                grid.v1Wind=grid.v2Wind;
            end
			

		end%main2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time step update of velocity fields from fvcompath.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		grid.unc1=grid.unc2;
		grid.vnc1=grid.vnc2;
		grid.wnc1=grid.wnc2;
		grid.elnc1=grid.elnc2;	
        if grid.diffusion
            grid.viscofhnc1=grid.viscofhnc2;
            grid.khnc1=grid.khnc2;
        end	
        if grid.wave
            grid.unc1Wave=grid.unc2Wave;
            grid.vnc1Wave=grid.vnc2Wave;
        end
        if grid.wind
            grid.unc1Wind=grid.unc2Wind;
            grid.vnc1Wind=grid.vnc2Wind;
        end

	end%mainfor



