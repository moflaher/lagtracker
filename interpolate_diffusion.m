function [lag,grid]=interpolate_diffusion(lag,grid,ifc) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: lag,grid
% Return: lag,grid
% Flags: 
%   ifc  = 0 if host is known to have correct host elements for current particle positions.    
%   ifc  = 1 if host should be updated (adds computational work particle positions.   
%   This is always zero because interpolate_diffusion is always currently called after interpolatev so the host position is known.   
%
%  Takes (xpt,ypt,sigpt) from lag and (viscofh,kh) from grid.
%  Obtains a linear interpolation of the provided field (viscofh,kh) at these points and there derivative in x,y for viscofh and z for kh.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%!==============================================================================!
%!  deTeRMine eLeMenT COnTaininG pOint (xp,yp,zp)                               !
%!==============================================================================!


  
  	if(ifc == 1)
		[allfound,lag]=newfindquick(lag,grid);
		if(allfound==0)
        		[lag,grid]=findfull(lag,grid,0);
     	end
  	end

%!===============================================================================!
%!  Find horizontal and vertical eddy viscosity at (xpt,ypt,sigpt)
%!===============================================================================!

    array=find(lag.indomain ~=0 & lag.inwater ~=0);    
    host=lag.host(array);
    sigpt=lag.sigpt(array);
	x0c=lag.xpt(array)-grid.xc(host);
	y0c=lag.ypt(array)-grid.yc(host);
    e=grid.nv(host,:);


  
    
    zf1=0*host;
    zf2=zf1;
    k1=ones(size(host));
    k2=k1;
    %    particle near surface
    a=sigpt>grid.zz(1);
        zf1(a)=1;
        zf2(a)=0;
        k1(a)=1;
        k2(a)=1;

    a=sigpt<grid.zz(grid.siglay);
        zf1(a)=1;
        zf2(a)=0;
        k1(a)=grid.siglay;
        k2(a)=grid.siglay;
        
    a=sigpt>grid.zz(grid.siglay)&sigpt<grid.zz(1);
        dzp=(1./grid.zz)*sigpt(a)';
        dzp(dzp<1)=10;
        [mtemp k1(a)] = min(dzp);
        k2(a)=k1(a)+1;		  
        zf1(a)=(sigpt(a)-grid.zz(k2(a)))./(grid.zz(k1(a))-grid.zz(k2(a)));
        zf2(a)=(grid.zz(k1(a))-sigpt(a))./(grid.zz(k1(a))-grid.zz(k2(a)));
    
        [m n]=size(grid.viscofhin);
        viscofh01 = grid.aw0(host,1).*grid.viscofhin(e(:,1)+m*(k1-1))+grid.aw0(host,2).*grid.viscofhin(e(:,2)+m*(k1-1))+grid.aw0(host,3).*grid.viscofhin(e(:,3)+m*(k1-1));
		viscofhx1 = grid.awx(host,1).*grid.viscofhin(e(:,1)+m*(k1-1))+grid.awx(host,2).*grid.viscofhin(e(:,2)+m*(k1-1))+grid.awx(host,3).*grid.viscofhin(e(:,3)+m*(k1-1));
		viscofhy1 = grid.awy(host,1).*grid.viscofhin(e(:,1)+m*(k1-1))+grid.awy(host,2).*grid.viscofhin(e(:,2)+m*(k1-1))+grid.awy(host,3).*grid.viscofhin(e(:,3)+m*(k1-1));
		viscofhe01 = viscofh01 + viscofhx1.*x0c + viscofhy1.*y0c;

        viscofh02 = grid.aw0(host,1).*grid.viscofhin(e(:,1)+m*(k2-1))+grid.aw0(host,2).*grid.viscofhin(e(:,2)+m*(k2-1))+grid.aw0(host,3).*grid.viscofhin(e(:,3)+m*(k2-1));
		viscofhx2 = grid.awx(host,1).*grid.viscofhin(e(:,1)+m*(k2-1))+grid.awx(host,2).*grid.viscofhin(e(:,2)+m*(k2-1))+grid.awx(host,3).*grid.viscofhin(e(:,3)+m*(k2-1));
		viscofhy2 = grid.awy(host,1).*grid.viscofhin(e(:,1)+m*(k2-1))+grid.awy(host,2).*grid.viscofhin(e(:,2)+m*(k2-1))+grid.awy(host,3).*grid.viscofhin(e(:,3)+m*(k2-1));
		viscofhe02 = viscofh02 + viscofhx2.*x0c + viscofhy2.*y0c;

        lag.viscofhp(array)= viscofhe01.*zf1 + viscofhe02.*zf2;
        lag.viscofhx(array)= viscofhx1.*zf1 + viscofhx2.*zf2;
        lag.viscofhy(array)= viscofhy1.*zf1 + viscofhy2.*zf2;



    a=sigpt==grid.z(1);
        zf1(a)=1;
        zf2(a)=0;
        k1(a)=1;
        k2(a)=1+1;
    a=sigpt==grid.z(grid.siglev);
        zf1(a)=0;
        zf2(a)=1;
        k1(a)=grid.siglev-1;
        k2(a)=grid.siglev;

    a=sigpt<grid.z(1) & sigpt>grid.z(grid.siglev);
        dzp=(1./grid.z)*sigpt(a)';
        dzp(dzp<1)=10;
        [mtemp k1(a)] = min(dzp);
        k2(a)=k1(a)+1;		  
        zf1(a)=(sigpt(a)-grid.z(k2(a)))./(grid.z(k1(a))-grid.z(k2(a)));
        zf2(a)=(grid.z(k1(a))-sigpt(a))./(grid.z(k1(a))-grid.z(k2(a)));



        [m n]=size(grid.khin);
        kh0 = grid.aw0(host,1).*grid.khin(e(:,1)+m*(k1-1))+grid.aw0(host,2).*grid.khin(e(:,2)+m*(k1-1))+grid.aw0(host,3).*grid.khin(e(:,3)+m*(k1-1));
		khx = grid.awx(host,1).*grid.khin(e(:,1)+m*(k1-1))+grid.awx(host,2).*grid.khin(e(:,2)+m*(k1-1))+grid.awx(host,3).*grid.khin(e(:,3)+m*(k1-1));
		khy = grid.awy(host,1).*grid.khin(e(:,1)+m*(k1-1))+grid.awy(host,2).*grid.khin(e(:,2)+m*(k1-1))+grid.awy(host,3).*grid.khin(e(:,3)+m*(k1-1));
		khe01 = kh0 + khx.*x0c + khy.*y0c;

        kh0 = grid.aw0(host,1).*grid.khin(e(:,1)+m*(k2-1))+grid.aw0(host,2).*grid.khin(e(:,2)+m*(k2-1))+grid.aw0(host,3).*grid.khin(e(:,3)+m*(k2-1));
		khx = grid.awx(host,1).*grid.khin(e(:,1)+m*(k2-1))+grid.awx(host,2).*grid.khin(e(:,2)+m*(k2-1))+grid.awx(host,3).*grid.khin(e(:,3)+m*(k2-1));
		khy = grid.awy(host,1).*grid.khin(e(:,1)+m*(k2-1))+grid.awy(host,2).*grid.khin(e(:,2)+m*(k2-1))+grid.awy(host,3).*grid.khin(e(:,3)+m*(k2-1));
		khe02 = kh0 + khx.*x0c + khy.*y0c;

        lag.khp(array)= khe01.*zf1 + khe02.*zf2;
  
        
        leveldiff=grid.hele(host).*grid.z(k2)-grid.hele(host).*grid.z(k1);
        lag.khz(array)= (khe02-khe01)./leveldiff;
        


    %Apply a noslip boundary at the bottom, setting all diffusion to zero if the particle is at the bottom
	lag.viscofhp(lag.sigpt==-1)=0; 
	lag.khp(lag.sigpt==-1)=0; 


end
    

