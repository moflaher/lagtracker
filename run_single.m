%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
% Return: 
%
%   Run code on single core and save particle postions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


start_times=12;

%single cpu 
	for starti=0:3:12
		starti
        clear savelag
        savelag=main(starti);
        savelag=struct(savelag)
        save(['savedir/kit4_kelp_20m_0.018/element_77567_s' num2str(starti) '.mat'],'savelag','-v7.3');
	end

