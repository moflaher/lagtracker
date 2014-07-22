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
           
        xlocl=200000;     
        xlocu=250000;
        ylocl=200000;
        ylocu=250000; 



        idx=find(grid.xc>= xlocl & grid.xc <=xlocu & grid.yc>= ylocl & grid.yc <=ylocu );
        foundelements=length(idx)
        array=1:100:foundelements;		

        set.xstart(:,1) = grid.xc(idx(array));
		set.ystart(:,1)  = grid.yc(idx(array));
		set.sigstart(:,1) = ones(length(idx(array)),1)*-0.0001;
   
elseif set.option==5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start particles at elements in mat-file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    if ~exist('cage_elements_sfm6_musq2_locations.mat')
       load('cage_elements_sfm6_musq2.mat')

        savex=[];
        savey=[];
        savesig=[];

        pp=5;
        xbox_size=200;
	    ybox_size=200;
	    num_parts=pp*100;
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
        save('cage_elements_sfm6_musq2_locations.mat','savex','savey','savesig')

      else
        display('Loading particle start location data.')
        load('cage_elements_sfm6_musq2_locations.mat')
      end


        set.xstart(:,1) = savex;
		set.ystart(:,1)  =  savey;
		set.sigstart(:,1) = savesig;
   
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start particles at elements centers for entire grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	set.xstart  = grid.xc;
	set.ystart = grid.yc;
	set.sigstart(1:grid.nele,1) = zeros(1:grid.nele,1);
	

end
