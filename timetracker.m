function [time]=timetracker(set)

time=timetrackerdef;

%if mod(loopsperhour,outtimesperhour) ~=0 then there will be no output
time.loopsperhour=set.internalstep;
time.outtimesperhour=set.outputstep;

ncfvcom=netcdf( set.fvcompath ,'nowrite');

time.starthour=set.start;
time.trackingtime=set.track;

%dti 	-> internal stepping
%instp 	-> stepping of input data
%dtout 	-> output stepping
time.instp=(ncfvcom{'time'}(2)-ncfvcom{'time'}(1));
time.dti=time.instp/time.loopsperhour;
time.dtout=time.instp/time.outtimesperhour;

time.finishhour=time.starthour+time.trackingtime;
time.itout=0;
time.iint=0;

time.i2=floor(time.instp/time.dti);
time.int2=floor(time.dtout/time.dti);

%for setting size of x,y,z,u,v,w arrays  
time.outt=((time.trackingtime)*time.outtimesperhour)+1;

close(ncfvcom);

