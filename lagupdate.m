function [lag,grid,time]=lagupdate(lag,grid,time)

%!==============================================================================|
%!  update particle positions, calculate scalar fields and particle velocities  |
%!==============================================================================|

	grid.u1=grid.unc1;
	grid.v1=grid.vnc1;
	grid.w1=grid.wnc1;
    grid.el1=grid.elnc1;
	

%!------------------------------------------------------------------------------|
%!   loop over the tracking period                                              |
%!------------------------------------------------------------------------------|

	for nh=time.starthour+1:time.finishhour%mainfor
		

%!--hourly reading of velocity fields in netcdf file
		[grid]=ncdread(grid,nh);
     
%!------------------------------------------------------------------------------|
%!   loop within one hour
%!------------------------------------------------------------------------------|
		for it=1:time.i2%main2

%			['On step: ' int2str(time.iint+1) '/' int2str(time.i2*time.trackingtime) ]

			time.iint=time.iint+1;

%!--time interpolation of physical fields
        		tmp2=it/time.i2;
        		tmp1=(it-time.i2)/-time.i2;

        		grid.u2=tmp1*grid.unc1+tmp2*grid.unc2;
                grid.v2=tmp1*grid.vnc1+tmp2*grid.vnc2;
        		grid.w2=tmp1*grid.wnc1+tmp2*grid.wnc2;
        		grid.el2=tmp1*grid.elnc1+tmp2*grid.elnc2;
                grid.dedtin=(grid.el2-grid.el1)/time.dti;
                
                thm1=mod(nh-1,grid.size)+1;
                th=mod(nh,grid.size)+1;
                lag.turbine_sigma=tmp1*grid.zz(grid.turbine_sigmas(thm1,:))+tmp2*grid.zz(grid.turbine_sigmas(th,:));
%!--particle is moving
			[lag grid time]=rungekutta(lag,grid,time);

%!--write data to file 
        		if(mod(time.iint,time.int2) == 0)%if3
%            			lag.up(:) =   (lag.xp(:)-lag.xplst(:))/time.dtout;
%            			lag.vp(:) =   (lag.yp(:)-lag.yplst(:))/time.dtout;
%            			lag.wp(:) =   ((lag.sigp(:)-lag.sigplst(:))/time.dtout).*(lag.hp+lag.ep);
                
                
%!--write particle position and velocities to file 
				time.itout=time.itout+1;

                array=lag.indomain==1;
                lag.x(array,time.itout)=lag.xp(array) + grid.vxmin;
           		lag.y(array,time.itout)=lag.yp(array) + grid.vymin;
				lag.sig(array,time.itout)=lag.sigp(array);
				lag.z(array,time.itout)=lag.sigp(array).*(lag.hp(array)+lag.ep(array))+lag.ep(array);
				lag.u(array,time.itout)=lag.up(array);
           		lag.v(array,time.itout)=lag.vp(array);
				lag.w(array,time.itout)=lag.wp(array);
                array=lag.indomain==0;
                lag.x(array,time.itout)=NaN;
           		lag.y(array,time.itout)=NaN;
				lag.sig(array,time.itout)=NaN;
				lag.z(array,time.itout)=NaN;
				lag.u(array,time.itout)=0;
           		lag.v(array,time.itout)=0;
				lag.w(array,time.itout)=0;
				lag.time(time.itout)=(time.iint*time.dti)+(time.starthour*time.instp);
				lag.xplst = lag.xp;
				lag.yplst = lag.yp;
				lag.sigplst = lag.sigp;
        		end%if3

%!--time step update of velocity fields
       		grid.u1=grid.u2;
			grid.v1=grid.v2;
			grid.w1=grid.w2;
			grid.el1=grid.el2;
			

		end%main2

%!--hourly update of physical fields
		grid.unc1(:)=grid.unc2(:);
		grid.vnc1(:)=grid.vnc2(:);
		grid.wnc1(:)=grid.wnc2(:);
		grid.elnc1(:)=grid.elnc2(:);		

	end%mainfor



