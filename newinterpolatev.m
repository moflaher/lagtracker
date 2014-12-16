function [lag,grid]=newinterpolatev(lag,grid,ifc) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: lag,grid
% Return: lag,grid
%
%  Takes (xpt,ypt,sigpt) from lag and (uin,vin,win) from grid.
%  Obtains a linear interpolation of the provided velocity field (uin,vin,win) at these points.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine element containing point (xp,yp,sigpt)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if(ifc == 1)
		[allfound,lag]=newfindquick(lag,grid);
		if(allfound==0)
        		[lag,grid]=findfull(lag,grid,0);
     	end
  	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine velocity at (xp,yp,sigpt)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    array=find(lag.indomain ==1 & lag.inwater ==1);
    
    host=lag.host(array);
    sigpt=lag.sigpt(array);
	x0c=lag.xpt(array)-grid.xc(host);
	y0c=lag.ypt(array)-grid.yc(host);
    e=grid.nbe(host,:);
    e(e==0)=grid.nele+1;
    
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
%         dzp=grid.zz*ones(size(sigpt(a)'))-ones(size(grid.zz))*sigpt(a)';
%         dzp(dzp<0)=1;

        if length(a)==1
            dzp=(1./grid.zz)*sigpt;    
        else
            dzp=(1./grid.zz)*sigpt(a)';
        end

        dzp(dzp<1)=10;
        [mtemp k1(a)] = min(dzp);
        k2(a)=k1(a)+1;		  
        zf1(a)=(sigpt(a)-grid.zz(k2(a)))./(grid.zz(k1(a))-grid.zz(k2(a)));
        zf2(a)=(grid.zz(k1(a))-sigpt(a))./(grid.zz(k1(a))-grid.zz(k2(a)));

    
    [m n]=size(grid.uin);
    dudx= grid.a1u(host,1).*grid.uin(host+m*(k1-1))+grid.a1u(host,2).*grid.uin(e(:,1)+m*(k1-1))+grid.a1u(host,3).*grid.uin(e(:,2)+m*(k1-1))+grid.a1u(host,4).*grid.uin(e(:,3)+m*(k1-1));
	dudy= grid.a2u(host,1).*grid.uin(host+m*(k1-1))+grid.a2u(host,2).*grid.uin(e(:,1)+m*(k1-1))+grid.a2u(host,3).*grid.uin(e(:,2)+m*(k1-1))+grid.a2u(host,4).*grid.uin(e(:,3)+m*(k1-1));
	ue01= grid.uin(host+m*(k1-1)) + dudx.*x0c + dudy.*y0c;
    
    dudx= grid.a1u(host,1).*grid.uin(host+m*(k2-1))+grid.a1u(host,2).*grid.uin(e(:,1)+m*(k2-1))+grid.a1u(host,3).*grid.uin(e(:,2)+m*(k2-1))+grid.a1u(host,4).*grid.uin(e(:,3)+m*(k2-1));
	dudy= grid.a2u(host,1).*grid.uin(host+m*(k2-1))+grid.a2u(host,2).*grid.uin(e(:,1)+m*(k2-1))+grid.a2u(host,3).*grid.uin(e(:,2)+m*(k2-1))+grid.a2u(host,4).*grid.uin(e(:,3)+m*(k2-1));
	ue02= grid.uin(host+m*(k2-1)) + dudx.*x0c + dudy.*y0c;

    lag.up(array)= ue01.*zf1 + ue02.*zf2;



    dvdx= grid.a1u(host,1).*grid.vin(host+m*(k1-1))+grid.a1u(host,2).*grid.vin(e(:,1)+m*(k1-1))+grid.a1u(host,3).*grid.vin(e(:,2)+m*(k1-1))+grid.a1u(host,4).*grid.vin(e(:,3)+m*(k1-1));
	dvdy= grid.a2u(host,1).*grid.vin(host+m*(k1-1))+grid.a2u(host,2).*grid.vin(e(:,1)+m*(k1-1))+grid.a2u(host,3).*grid.vin(e(:,2)+m*(k1-1))+grid.a2u(host,4).*grid.vin(e(:,3)+m*(k1-1));
	ve01= grid.vin(host+m*(k1-1)) + dvdx.*x0c + dvdy.*y0c;
    
    dvdx= grid.a1u(host,1).*grid.vin(host+m*(k2-1))+grid.a1u(host,2).*grid.vin(e(:,1)+m*(k2-1))+grid.a1u(host,3).*grid.vin(e(:,2)+m*(k2-1))+grid.a1u(host,4).*grid.vin(e(:,3)+m*(k2-1));
	dvdy= grid.a2u(host,1).*grid.vin(host+m*(k2-1))+grid.a2u(host,2).*grid.vin(e(:,1)+m*(k2-1))+grid.a2u(host,3).*grid.vin(e(:,2)+m*(k2-1))+grid.a2u(host,4).*grid.vin(e(:,3)+m*(k2-1));
	ve02= grid.vin(host+m*(k2-1)) + dvdx.*x0c + dvdy.*y0c;
    
    lag.vp(array)= ve01.*zf1 + ve02.*zf2;

    %a=sigpt==grid.z(1);
    %    zf1(a)=1;
    %    zf2(a)=0;
    %    k1(a)=1;
    %    k2(a)=1;
    %a=sigpt==grid.z(grid.siglev);
    %    zf1(a)=0;
    %    zf2(a)=1;
    %    k1(a)=grid.siglev;
    %    k2(a)=grid.siglev;

    %a=sigpt<grid.z(1) & sigpt>grid.z(grid.siglev);
%    %    dzp=grid.z*ones(size(sigpt(a)'))-ones(size(grid.z))*sigpt(a)';
%     %   dzp(dzp<0)=1;
      %  dzp=(1./grid.z)*sigpt(a)';
      %  dzp(dzp<1)=10;
      %  [mtemp k1(a)] = min(dzp);
      %  k2(a)=k1(a)+1;		  
      %  zf1(a)=(sigpt(a)-grid.z(k2(a)))./(grid.z(k1(a))-grid.z(k2(a)));
      %  zf2(a)=(grid.z(k1(a))-sigpt(a))./(grid.z(k1(a))-grid.z(k2(a)));

    [m n]=size(grid.win);

    dwdx= grid.a1u(host,1).*grid.win(host+m*(k1-1))+grid.a1u(host,2).*grid.win(e(:,1)+m*(k1-1))+grid.a1u(host,3).*grid.win(e(:,2)+m*(k1-1))+grid.a1u(host,4).*grid.win(e(:,3)+m*(k1-1));
	dwdy= grid.a2u(host,1).*grid.win(host+m*(k1-1))+grid.a2u(host,2).*grid.win(e(:,1)+m*(k1-1))+grid.a2u(host,3).*grid.win(e(:,2)+m*(k1-1))+grid.a2u(host,4).*grid.win(e(:,3)+m*(k1-1));
	we01= grid.win(host+m*(k1-1)) + dwdx.*x0c + dwdy.*y0c;
    
    dwdx= grid.a1u(host,1).*grid.win(host+m*(k2-1))+grid.a1u(host,2).*grid.win(e(:,1)+m*(k2-1))+grid.a1u(host,3).*grid.win(e(:,2)+m*(k2-1))+grid.a1u(host,4).*grid.win(e(:,3)+m*(k2-1));
	dwdy= grid.a2u(host,1).*grid.win(host+m*(k2-1))+grid.a2u(host,2).*grid.win(e(:,1)+m*(k2-1))+grid.a2u(host,3).*grid.win(e(:,2)+m*(k2-1))+grid.a2u(host,4).*grid.win(e(:,3)+m*(k2-1));
	we02= grid.win(host+m*(k2-1)) + dwdx.*x0c + dwdy.*y0c;
    
    lag.wp(array)= we01.*zf1 + we02.*zf2;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply a noslip boundary at the bottom, setting all velocities to zero if the particle is at the bottom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	lag.up(lag.sigpt==-1)=0; 
	lag.vp(lag.sigpt==-1)=0; 
	lag.wp(lag.sigpt==-1)=0;

end
    

