%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
% Return: 
%
%   Run code on single core and save particle postions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


start_times=12;

%single cpu 
	for starti=0:3:start_times
		starti
		savelag=main(starti);
        savelag=struct(savelag)
        save(['savedir/25_part_no_cage_in60min_time1min_out10min/no_cages_25_part_sfm6_musq2_' num2str(starti) '.mat'],'savelag');
	end

