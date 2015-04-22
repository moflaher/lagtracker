function [lag,grid]=findfull(lag,grid,pter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: lag,grid
% Return: lag,grid
% Flags: pter
%   pter==1 sets initial particle positions
%   pter==0 all other times
% 
%  Determine Which Element A List of Particles Reside in By Searching Nearest Element to Progressively Further Elements.  
%  Updates host component of lagrangian Particle. 
%  Updates array "ifound" flagging whether the host has been found.	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



	if pter==0
		xp=lag.xpt;
		yp=lag.ypt; 
	else
		xp=lag.xp;
		yp=lag.yp;
    end

    	array=find(lag.ifound ==0 & lag.indomain~=0 & lag.inwater~=0);
        thost=nan(length(xp),1);
        idxa=find(grid.vx>= min(min(xp)) &grid.vx<= max(max(xp)) &grid.vy>= min(min(yp)) &grid.vy<= max(max(yp)));
        idxa=ismember(grid.nv,idxa);
        uvidx2=1:length(idxa);
        idxa=uvidx2(logical(sum(idxa')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two methods for finding host.
%    First is based on finding host of each particle (original).
%    Second is based on finding which particles are in each host.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if true%length(array)<=length(idxa)
	    for jj=1:length(array)%mainfor
		    i=array(jj);
			xlag=xp(i);
		    ylag=yp(i);
		    radlist(1:grid.nele,1)=(grid.xc(1:grid.nele)-xlag).^2 + (grid.yc(1:grid.nele)-ylag).^2;
			%radlist=(grid.xc-xlag).^2 + (grid.yc-ylag).^2;
		    %check closest triangle
		    [minrad mini]=min(radlist);
		    xtri=grid.vx(grid.nv(mini,1:3))';
		    ytri=grid.vy(grid.nv(mini,1:3))'; 
		    if(multitrianglechecker(xtri,ytri,xlag,ylag)==1)%if2
			    lag.ifound(i)=1; 
			    lag.host(i)=mini;	
		    else%if2
			    %tri=find(multitrianglechecker(grid.vx(grid.nv(:,1:3)),grid.vy(grid.nv(:,1:3)),xlag,ylag));
				tri=find(multitrianglechecker(grid.vx(grid.nv),grid.vy(grid.nv),xlag,ylag));
			    if ~isempty(tri)
		                    lag.ifound(i)=1; 
		                    lag.host(i)=tri(1);
%		                    if length(tri)>1 
%		                        tri
%		                    end
                %else
                %'no find' 
                %i
				end                
            end%if2

	    end%mainfor
        sum(lag.ifound)/length(lag.ifound)
	    lag.indomain=0*lag.indomain;
	    lag.indomain(lag.ifound == 1)=1;
	    lag.inwater=0*lag.inwater;
	    lag.inwater(lag.ifound == 1)=1;
    else
        thost=nan(length(xp),1);
        idxa=find(grid.vx>= min(min(xp)) &grid.vx<= max(max(xp)) &grid.vy>= min(min(yp)) &grid.vy<= max(max(yp)));
        idxa=ismember(grid.nv,idxa);
        uvidx2=1:length(idxa);
        idxa=uvidx2(logical(sum(idxa')));

        for i=1:length(idxa)
            thost(inpolygon(xp,yp,grid.vx(grid.nv(idxa(i),1:3)),grid.vy(grid.nv(idxa(i),1:3))))=idxa(i);   
        end    

        lag.ifound(~isnan(thost))=1;
        lag.host(~isnan(thost))=thost(~isnan(thost));

	    lag.indomain=0*lag.indomain;
	    lag.indomain(lag.ifound == 1)=1;
	    lag.inwater=0*lag.inwater;
	    lag.inwater(lag.ifound == 1)=1;

    end

end
