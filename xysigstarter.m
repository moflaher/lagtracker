function [grid,time,set]=xysigstarter(grid,time,set)


if set.option==1
	% Put a random box of n particles around each turbine 
	xbox_size=20;
	ybox_size=20;
	num_parts=500;
	for i=1:grid.nturbines
		xs=(rand(num_parts,1)-0.5)*xbox_size;
		ys=(rand(num_parts,1)-0.5)*ybox_size;
		sigs= -1 + rand(num_parts,1)*(1+grid.zz(grid.turbine_sigmas(i,1)+1));
		set.xstart(1+(i-1)*num_parts:i*num_parts,1) = grid.xc(grid.turbines(i))+xs;
		set.ystart(1+(i-1)*num_parts:i*num_parts,1)  = grid.yc(grid.turbines(i))+ys;
		set.sigstart(1+(i-1)*num_parts:i*num_parts,1) = sigs;
    	end
else
	%starts particle in the center of each element at a sigma lvl	
	set.xstart  = grid.xc;
	set.ystart = grid.yc;
	set.sigstart(1:grid.nele,1) = 0;
	

end
