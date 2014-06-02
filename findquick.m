function [allfound,lag,grid]=findquick(lag,grid)%,pter)
%!==============================================================================|
%!  Determine Which Element A List of Particles Reside in By Searching          |
%!  Neighboring Elements.  Updates host component of lagrangian Particle        |  
%!  Type and updates logical array "ELEM_found" flagging whether the host       |
%!  Has been found							         |
%!==============================================================================|

	allfound=0;
    lag.ifound=0*lag.ifound;
    lag.sbound=0*lag.sbound;
    
	for i=1:lag.npts		%mainfor
		if(lag.indomain(i)~=0) 	%mainif
			%pter=0;
			%if pter==0
			%	xlag=lag.xpt(i);
			%	ylag=lag.ypt(i); 
			%else
			%	xlag=lag.xp(i);
			%	ylag=lag.yp(i);
			%end
			xlag=lag.xpt(i);
			ylag=lag.ypt(i); 
			ilast=lag.host(i);
			xlast=grid.vx(grid.nv(ilast,1:3));
			ylast=grid.vy(grid.nv(ilast,1:3));
			%particle remains in element
			if(trianglechecker(xlast,ylast,xlag,ylag)==1)%if2   
				lag.ifound(i)=1;
			%check neighbors     
		    else                                                   
				for j=1:3%for2
					ncheck=grid.nv(ilast,j);
					for k=1:grid.ntve(grid.nv(ilast,j))%for3
						iney=grid.nbve(ncheck,k); 
						xney=grid.vx(grid.nv(iney,1:3));
						yney=grid.vy(grid.nv(iney,1:3));
						if(trianglechecker(xney,yney,xlag,ylag)==1)
							lag.ifound(i)=1; 
							lag.host(i)=iney;
							if(grid.isonb(grid.nv(iney,1))==1 || grid.isonb(grid.nv(iney,2))==1 || grid.isonb(grid.nv(iney,3))==1)
								lag.sbound(i)=1;
							end
							break;  
						end
					end%for3
					if (lag.ifound(i)==1)
					break;
					end
				end%for2
			end%if2
		end%mainif
	end%mainfor
	  
	if(sum(lag.ifound) == sum(lag.indomain)) 
		allfound=1;
	end


