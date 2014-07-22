function [allfound,lag]=newfindquick(lag,grid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: lag,grid
% Return: lag,allfound
% 
%  Determine Which Element A List of Particles Reside in By Searching Neighboring Elements.  
%  Updates host component of lagrangian Particle. 
%  Updates array "ifound" flagging whether the host has been found.	
%  allfound shows if the locations of all particles are known.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    lag.ifound=0*lag.ifound;
    lag.sbound=0*lag.sbound;
    array=find(lag.indomain~=0);
	
	xlag=lag.xpt(array); 
	ylag=lag.ypt(array); 
	ilast=lag.host(array);
    xlasttri=grid.vx(grid.nv(ilast,1:3));
	ylasttri=grid.vy(grid.nv(ilast,1:3));
    [aa bb]=size(xlasttri);
    if bb==1
        xlasttri=xlasttri';
        ylasttri=ylasttri';
    end
	lag.ifound(array)=multitrianglechecker(xlasttri,ylasttri,xlag,ylag);%if2   
	
    [m n]=size(grid.closest);
    for nbi=1:n
        array=find(lag.indomain~=0 & lag.ifound==0);
        if isempty(array), break; end;
        ibe=grid.closest(lag.host(array),nbi);
        a=find(ibe>0);
        a2=array(a);
        ibe=ibe(a);
        xlag=lag.xpt(a2); 
        ylag=lag.ypt(a2);        
        xbe=grid.vx(grid.nv(ibe,1:3));
        ybe=grid.vy(grid.nv(ibe,1:3)); 
        if length(a)==1
            xbe=xbe';
            ybe=ybe';
        end
        lag.ifound(a2)=multitrianglechecker(xbe,ybe,xlag,ylag);%if2
        a3=lag.ifound(a2)==1;
        lag.host(a2(a3))=ibe(a3);
    end
    allfound=(sum(lag.ifound) == sum(lag.indomain));
end


