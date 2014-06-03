function [time]=timetracker(set)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: set
% Return: time
%
%   Sets up time settings.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


time=timetrackerdef;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if mod(loopsperhour,outtimesperhour) ~=0 then there will be no output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time.loopsperhour=set.internalstep;
time.outtimesperhour=set.outputstep;

ncid=netcdf.open(set.fvcompath,'NC_NOWRITE');

time.starthour=set.start;
time.trackingtime=set.track;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dti 	-> internal stepping
%instp 	-> stepping of input data
%dtout 	-> output stepping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ttime = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time')));
[Y, M, D, H, MN, S]=datevec(ttime(3)-ttime(2));
tsec=(H*3600)+(MN*60)+S;

time.instp=tsec;

time.dti=time.instp/time.loopsperhour;
time.dtout=time.instp/time.outtimesperhour;

time.finishhour=time.starthour+time.trackingtime;
time.itout=0;
time.iint=0;

time.i2=floor(time.instp/time.dti);
time.int2=floor(time.dtout/time.dti);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for setting size of x,y,z,u,v,w arrays  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time.outt=((time.trackingtime)*time.outtimesperhour)+1;

netcdf.close(ncid);
