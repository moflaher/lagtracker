function output=main_turbine_box(seti,swim_speed)

%grabs settings
	[set]=settings(seti,swim_speed);

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

    if seti==1, save gridfile grid; end
	    hosts=lag.host;
	    distance_travelled=0*hosts; 
	    for i=1:length(distance_travelled)
		array=find(lag.sig(i,:)==-1);
		if isempty(array)
		    distance_travelled(i)=NaN;
		    hosts(i)=NaN;
		else
		    bot_time=min(array);
		    distance_travelled(i)=sqrt((lag.x(i,bot_time)-lag.x(i,1)).^2+(lag.y(i,bot_time)-lag.y(i,1)).^2);
		end
	    end
	    
	    output=hosts;
	    output(lag.npts+1:2*lag.npts)=distance_travelled;
	
end
    



