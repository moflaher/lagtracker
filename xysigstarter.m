function [grid,time,set]=xysigstarter(grid,time,set)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: grid,time,set
% Return: grid,time,set
%
%   Sets up particle initial locations.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if set.option==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put a random box of n particles around each turbine 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
elseif set.option==2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put a random box of n particles around a element number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	xbox_size=5000;
	ybox_size=5000;
	num_parts=10;
	
		xs=(rand(num_parts,1)-0.5)*xbox_size;
		ys=(rand(num_parts,1)-0.5)*ybox_size;
        %loc=206894;
        loc=45177;
		set.xstart(:,1) = grid.xc(loc)+xs;
		set.ystart(:,1)  = grid.yc(loc)+ys;
		set.sigstart(:,1) = ones(num_parts,1)*-0.0001;

elseif set.option==3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put a random box of n particles around a location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	xbox_size=5000;
	ybox_size=5000;
	num_parts=1000;
	
		xs=(rand(num_parts,1)-0.5)*xbox_size;
		ys=(rand(num_parts,1)-0.5)*ybox_size;
         
        %kitimat         
        %xloc=242600;
        %yloc=310000;
        %gilislandleft         
        %xloc=193800;
        %yloc=225500;
        %gilislandright         
        %xloc=213000;
        %yloc=225500;
        %gilislandsouth         
        xloc=190000;
        yloc=190000;
		
        set.xstart(:,1) = xloc+xs;
		set.ystart(:,1)  = yloc+ys;
		set.sigstart(:,1) = zeros(num_parts,1);
   
elseif set.option==4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start particles at elements centers found in specified box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
		
        %xlocl=433000;     
        %xlocu=436000;
        %ylocl=729000;
        %ylocu=736000;       
           

        %model xy gets converted to lagtracker xy
        xlocl=-405525.84375;     
        xlocu=138350.4375;
        ylocl=-519105.15625;
        ylocu=-151528.140625; 


        idx=find(grid.xc>= xlocl-grid.vxmin & grid.xc <=xlocu-grid.vxmin & grid.yc>= ylocl-grid.vymin & grid.yc <=ylocu-grid.vymin );
        foundelements=length(idx)
        array=1:1:foundelements;		

        set.xstart(:,1) = grid.xc(idx(array));
		set.ystart(:,1)  = grid.yc(idx(array));
		set.zstart(:,1) = zeros(length(idx(array)),1);
   
elseif set.option==5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start particles at elements in mat-file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    elementfile='element_starts/cage_elements_sfm6_musq2.mat'
    locationfile='element_starts/cage_elements_sfm6_musq2_locations_10pp.mat'

    if ~exist(locationfile)
       load(elementfile)

        savex=[];
        savey=[];
        savesig=[];

        pp=10;
        xbox_size=300;
	    ybox_size=300;
	    num_parts=pp*1000;
        for i=1:length(cage_elements)

		    xs=(rand(num_parts,1)-0.5)*xbox_size;
		    ys=(rand(num_parts,1)-0.5)*ybox_size;
		    sigs= -1.*rand(num_parts,1);
            sigs=sigs(sigs>=-1 & sigs<=0);
	        xs=grid.xc(cage_elements(i))+xs;
		    ys=grid.yc(cage_elements(i))+ys;
		    host=inpolygon(xs,ys,grid.vx(grid.nv(cage_elements(i),:)),grid.vy(grid.nv(cage_elements(i),:)));
            
            xs=xs(host);
            ys=ys(host);
            sigs=sigs(host);
            savex=[savex; xs(1:pp)];
            savey=[savey; ys(1:pp)]; 
            savesig=[savesig; sigs(1:pp)];         

        end  	
        save(locationfile,'savex','savey')

      else
        display('Loading particle start location data.')
        load(locationfile)
      end


        set.xstart(:,1) = savex;
		set.ystart(:,1)  =  savey;
		set.zstart(:,1) = zeros(length(savex),1);
elseif set.option==6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start particles in box around elements in mat-file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    elementfile='element_starts/kit4_kelp_tight2_kelpfield_3elements.mat'
    locationfile='element_starts/kit4_kelp_tight2_kelpfield_3elements_200x200_1000pp.mat'

    if ~exist(locationfile)
       load(elementfile)

        savex=[];
        savey=[];

        pp=1000;
        xbox_size=200;
	    ybox_size=200;
        for i=1:length(elements)

		    xs=(rand(pp,1)-0.5)*xbox_size;
		    ys=(rand(pp,1)-0.5)*ybox_size;
		
	        xs=grid.xc(elements(i))+xs;
		    ys=grid.yc(elements(i))+ys;

            savex=[savex; xs];
            savey=[savey; ys]; 
  

        end  	
        save(locationfile,'savex','savey')

      else
        display('Loading particle start location data.')
        load(locationfile)
      end


        set.xstart(:,1) = savex;
		set.ystart(:,1)  =  savey;
		set.zstart(:,1) = zeros(length(savex),1);
   
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start particles at elements centers for entire grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	set.xstart  = grid.xc;
	set.ystart = grid.yc;
	set.zstart(1:grid.nele,1) = zeros(1:grid.nele,1);
	

end
