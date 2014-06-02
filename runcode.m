

%add mexcdf paths (if not already in your setpath)
addpath /array/home/user/static/mexcdf/netcdf_toolbox/netcdf/ncutility
addpath /array/home/user/static/mexcdf/netcdf_toolbox/netcdf/nctype
addpath /array/home/user/static/mexcdf/netcdf_toolbox/netcdf
addpath /array/home/user/static/mexcdf/mexnc




%runcode options
start_times=96;
number_of_elements=8629;
number_of_turbines=3;

%initialize storage array
turbineintersect=zeros(start_times,number_of_elements,number_of_turbines);


%single cpu 
	for starti=1:start_times
		starti
		turbineintersect(starti,:,:)=main(starti);
	end


%parallel cpu
%matlabpool open 7
%	parfor starti=1:start_times
%		starti
%		turbineintersect(starti,:,:)=main(starti);
%	end
%matlabpool close


save turbineintersect.mat turbineintersect









