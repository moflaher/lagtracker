
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
% Return: 
% 
%  Code to run main_turbine_box should be updated.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%number of particle released in a box
numparts=500;
%number of turbines
numturbs=3;

totalparts=numturbs*numparts;

%single cpu
for swimi=1:10
	%sink speed of a particle (sigma)
    	swim_speed(swimi)=-0.0002*swimi
	ns=120;
	output=zeros(ns,2*totalparts);
	
	for starti=1:ns
		starti
		output(starti,:)=main_turbine_box(starti,swim_speed(swimi));
	end
	
	strg=['output' int2str(swimi) ]
	save(strg, 'output')
end


%parallel cpu
%matlabpool open 7
%	parfor swimi=1:10
%		%sink speed of a particle (sigma)
%	    	swim_speed(swimi)=-0.0002*swimi
%		ns=120;
%		output=zeros(ns,2*totalparts);
%	
%		for starti=1:ns
%			starti
%			output(starti,:)=main_turbine_box(starti,swim_speed(swimi));
%		end
%	
%		strg=['output' int2str(swimi) ]
%		save(strg, 'output')
%	end
%matlabpool close


save particle_data.mat swim_speed average_distance all_particles_in_triangle









