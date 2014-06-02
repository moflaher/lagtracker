function [lag,grid]=newinterpolateelh(lag,grid,ifc) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: lag,grid
% Return: lag,grid
% Flags: 
%   ifc  = 0 if host is known to have correct host elements for current particle positions.    
%   ifc  = 1 if host should be updated (adds computational work particle positions.   
%   This is always zero because newinterpolateelh is always currently called after interpolatev so the host position is known.   
%
%  Takes (xpt,ypt) from lag and (hin,ein,dedtin) from grid.
%  Obtains a linear interpolation of the provided fields (hin,ein,dedtin) at these points.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine element containing point (xp,yp)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if(ifc == 1)
		[allfound,lag]=newfindquick(lag,grid);
		if(allfound==0)
        		[lag,grid]=findfull(lag,grid,0);
     	end
  	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linearly interpolate free surface height and bathymetry   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    array=find(lag.indomain ~=0 & lag.inwater ~=0);
    
    host=lag.host(array);
	x0c=lag.xpt(array)-grid.xc(host);
	y0c=lag.ypt(array)-grid.yc(host);
    n=grid.nv(host,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear interpolation of bathymetry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			h0 = grid.aw0(host,1).*grid.hin(n(:,1))+grid.aw0(host,2).*grid.hin(n(:,2))+grid.aw0(host,3).*grid.hin(n(:,3));
			hx = grid.awx(host,1).*grid.hin(n(:,1))+grid.awx(host,2).*grid.hin(n(:,2))+grid.awx(host,3).*grid.hin(n(:,3));
			hy = grid.awy(host,1).*grid.hin(n(:,1))+grid.awy(host,2).*grid.hin(n(:,2))+grid.awy(host,3).*grid.hin(n(:,3));
			lag.hp(array) = h0 + hx.*x0c + hy.*y0c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear interpolation of free surface height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			e0 = grid.aw0(host,1).*grid.ein(n(:,1))+grid.aw0(host,2).*grid.ein(n(:,2))+grid.aw0(host,3).*grid.ein(n(:,3));
			ex = grid.awx(host,1).*grid.ein(n(:,1))+grid.awx(host,2).*grid.ein(n(:,2))+grid.awx(host,3).*grid.ein(n(:,3));
			ey = grid.awy(host,1).*grid.ein(n(:,1))+grid.awy(host,2).*grid.ein(n(:,2))+grid.awy(host,3).*grid.ein(n(:,3));
			lag.ep(array) = e0 + ex.*x0c + ey.*y0c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear interpolation of time derivative of free surface height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            e0 = grid.aw0(host,1).*grid.dedtin(n(:,1))+grid.aw0(host,2).*grid.dedtin(n(:,2))+grid.aw0(host,3).*grid.dedtin(n(:,3));
			ex = grid.awx(host,1).*grid.dedtin(n(:,1))+grid.awx(host,2).*grid.dedtin(n(:,2))+grid.awx(host,3).*grid.dedtin(n(:,3));
			ey = grid.awy(host,1).*grid.dedtin(n(:,1))+grid.awy(host,2).*grid.dedtin(n(:,2))+grid.awy(host,3).*grid.dedtin(n(:,3));
			lag.dedtp(array) = e0 + ex.*x0c + ey.*y0c;

end




