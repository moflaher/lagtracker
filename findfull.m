function [lag,grid]=findfull(lag,grid,pter)

%!==============================================================================|
%!  Find Home Element For Points (X,Y)                                          |
%!  Search Nearest Element to Progressively Further Elements. Updates Lagrangian|  
%!  component "host" and marks Lagrangian component "ifound" with 1 if          |
%!  found.  Returns logical variable "all_found" if all lagrangian variables    |
%!  have a known host element.  The host element may have been found prior to   |
%!  entry in this routine.                                                      |
%!==============================================================================|

  

	if pter==0
		xp=lag.xpt;
		yp=lag.ypt; 
	else
		xp=lag.xp;
		yp=lag.yp;
    	end

    	array=find(lag.ifound ==0 & lag.indomain~=0);
	
	for jj=1:length(array)%mainfor
		i=array(jj);
        	xlag=xp(i);
		ylag=yp(i);
		radlist(1:grid.nele,1)=(grid.xc(1:grid.nele)-xlag).^2 + (grid.yc(1:grid.nele)-ylag).^2;
		%check closest triangle
		[minrad mini]=min(radlist);
		xtri=grid.vx(grid.nv(mini,1:3));
		ytri=grid.vy(grid.nv(mini,1:3)); 
		if(trianglechecker(xtri,ytri,xlag,ylag)==1)%if2
			lag.ifound(i)=1; 
			lag.host(i)=mini;	
		else%if2
			tri=newtrianglechecker(grid.vx(grid.nv(:,1:3)),grid.vy(grid.nv(:,1:3)),xlag,ylag);
			if ~isempty(tri)
		                lag.ifound(i)=1; 
		                lag.host(i)=tri(1);
		                if length(tri)>1 
		                    tri
		                end
                    	end                
           	 end%if2
	end%mainfor

	lag.indomain=0*lag.indomain;
	lag.indomain(lag.ifound == 1)=1;
	lag.inwater=0*lag.inwater;
	lag.inwater(lag.ifound == 1 & lag.sbound == 0)=1;
end
