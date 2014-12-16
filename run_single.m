%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
% Return: 
%
%   Run code on single core and save particle postions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


start_times=12;

%single cpu 
	for starti=0
		starti
        clear savelag
        savelag=main(starti);
        savelag=struct(savelag)
        save(['savedir/sfm6_musq2_no_cages/cage_elements_s' num2str(starti) '_sink_1mm.mat'],'savelag','-v7.3');
	end

