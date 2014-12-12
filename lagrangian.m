function [lag]=lagrangian(time,set)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: time,set
% Return: lag
%
%   Initializes arrays for lag.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



n=length(set.xstart);


lag=lagrangiandef;

%%%%%%% some of these are outdated and need to be updated.
%itag       !!label for the particle 
%host       !!element containing particle
%found      !!host element is found       
%indomain   !!particle is in the domain   
%sbound     !!host element has a solid boundary node   
%xp         !!x position of particle 
%yp         !!y position of particle 
%sigp       !!sigmaposition of particle 
%xpt        !!x absolute position of particle
%ypt        !!y absolute position of particle 
%sigpt      !!z absolute position of particle 
%up         !!u velocity of particle 
%vp         !!v velocity of particle 
%wp         !!omega velocity of particle
%hp         !!bathymetry at particle position
%ep         !!free surface height at particle 
%w_swimspeed !!swim_speed in z coords
%omega_swimspeed !!swim_speed in sigma coords


lag.npts=n;
lag.time=zeros(time.outt,1);
lag.host=zeros(n,1);
lag.ifound=zeros(n,1);
lag.indomain=zeros(n,1);
lag.sbound=zeros(n,1);
lag.inwater=zeros(n,1);
lag.chix=zeros(n,4);
lag.chiy=zeros(n,4);
lag.chiz=zeros(n,4);
lag.xp=zeros(n,1);
lag.yp=zeros(n,1);
lag.zp=zeros(n,1);
lag.sigp=zeros(n,1);
lag.xpt=zeros(n,1);
lag.ypt=zeros(n,1);
lag.zpt=zeros(n,1);
lag.sigpt=zeros(n,1);
lag.up=zeros(n,1);
lag.vp=zeros(n,1);
lag.wp=zeros(n,1);
lag.w_swimspeed=zeros(n,1);
lag.input_swimspeed=set.swimspeed;
lag.hp=zeros(n,1);
lag.ep=zeros(n,1);
lag.dedtp=zeros(n,1);

lag.x=zeros(n,time.outt);
lag.y=zeros(n,time.outt);
lag.sig=zeros(n,time.outt);
lag.z=zeros(n,time.outt);
lag.u=zeros(n,time.outt);
lag.v=zeros(n,time.outt);
lag.w=zeros(n,time.outt);
lag.h=zeros(n,time.outt);
lag.check_turbine_intersect=set.check_turbines;

lag.xpt=set.xstart;
lag.ypt=set.ystart;
lag.zpt=set.zstart;
%lag.sigpt=set.sigstart;




if set.diffusion
    lag.savediffusionfudgefactor=set.diffusionfudgefactor;
    lag.viscofhp=zeros(n,1);
    lag.khp=zeros(n,1);
	if set.randomstate
		lag.saverandomstate=rng;
		lag.wiener=sqrt(time.dti)*randn(4*(((time.trackingtime)*time.loopsperhour)+1),1);	
	else
		lag.saverandomstate=set.randomstate;
		rng(set.randomstate);
		lag.wiener=sqrt(time.dti)*randn(4*(((time.trackingtime)*time.loopsperhour)+1),1);
	end
	
end







