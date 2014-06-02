function turbineintersect=main(seti)

%grabs settings
	[set]=settings(seti,0);

%sets up time variables
	[time]=timetracker(set);

%sets grid details up
	[grid]=gridsetup(time,set);
	[grid]=triangleedgegrid(grid);


%sets up particle start locations
	if set.starter==1
		[grid,time,set]=xysigstarter(grid,time,set);
	end

%set up stuff
	[lag]=lagrangian(time,set);
	[lag,grid,time]=setlag(lag,grid,time);
    	[lag]=check_turbine_intersections(lag,grid);

%move the particles
	[lag,grid,time]=lagupdate(lag,grid,time);
	
    	turbineintersect=lag.turbine_intersects;
    
   	if seti==1, save gridfile grid; end

	


end
    


