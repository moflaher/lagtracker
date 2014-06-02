function [set]=settings(seti,swim_speed)

set=settingsdef;

set.fvcompath=['/home/user/Desktop/lagtracker_clean/upper1_0001.nc'];


%time settings
set.start=0;
set.track=4820;
set.internalstep=1;
set.outputstep=1;


%lagrangian settings



%grid settings
set.numberoftimesteps=96;


%xysigstart settings
set.starter=1;
set.option=3;
set.start=floor(seti-1);
set.swimspeed=swim_speed;
set.check_turbines=true;














