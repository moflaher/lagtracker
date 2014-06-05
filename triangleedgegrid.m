function [grid]=triangleedgegrid(grid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: grid
% Return: grid
%
%   Calculates information about grid needed to increase particle tracking speed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!==============================================================================!
%!     variable list:							         !
%!  vx(m)    :: vx(i) = x-coordinate of node i (input from mesh)	         !
%!  vy(m)    :: vy(i) = y-coordinate of node i (input from mesh)                !
%!  nv(n,3)  :: nv(i:1-3) = 3 node numbers of element i                         !
%!  xc(n)    :: xc(i) = x-coordinate of element i (calculated from vx)          !
%!  yc(n)    :: yc(i) = y-coordinate of element i (calculated from vy)          !
%!                                                                              !
%!  nbe(n,3) :: nbe(i,1->3) = element index of 1->3 neighbors of element i      !
%!  isbce(n) :: flag if element is on the boundary, see below for values        !
%!  isonb(m) :: flag is node is on the boundary, see below for values           !
%!==============================================================================!
%!     classification of the triangles nodes, and edges                         !
%!                                                                              !
%!     isonb(i)=0:  node in the interior computational domain                   !
%!     isonb(i)=1:  node on the solid boundary                                  !
%!     isonb(i)=2:  node on the open boundary                                   !
%!                                                                              !
%!     isbce(i)=0:  element in the interior computational domain                !
%!     isbce(i)=1:  element on the solid boundary                               !
%!     isbce(i)=2:  element on the open boundary                                !
%!     isbce(i)=3:  element with 2 solid boundary edges                         !
%!                                                                              !
%!==============================================================================!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Set up mesh (horizontal coordinates)                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculate global minimums and maximums                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	grid.vxmin = min(grid.vx(1:grid.node));
	grid.vxmax = max(grid.vx(1:grid.node));
	grid.vymin = min(grid.vy(1:grid.node));
	grid.vymax = max(grid.vy(1:grid.node));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Shift grid to upper right cartesian                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	grid.vx(:) = grid.vx(:) - grid.vxmin;
	grid.vy(:) = grid.vy(:) - grid.vymin;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculate global element center grid coordinates                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		grid.xc  = (grid.vx(grid.nv(:,1)) + grid.vx(grid.nv(:,2)) + grid.vx(grid.nv(:,3)))/3;
		grid.yc  = (grid.vy(grid.nv(:,1)) + grid.vy(grid.nv(:,2)) + grid.vy(grid.nv(:,3)))/3;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   If element on boundary set isbce(i)=1 and isonb(j)=1 for boundary nodes                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	grid.isbce(grid.nbe(:,1)==0| grid.nbe(:,2)==0 | grid.nbe(:,3)==0) = 1;
	grid.isonb(grid.nv(grid.nbe(:,1) == 0,2)) = 1;
	grid.isonb(grid.nv(grid.nbe(:,1) == 0,3)) = 1;
	grid.isonb(grid.nv(grid.nbe(:,2) == 0,1)) = 1;
	grid.isonb(grid.nv(grid.nbe(:,2) == 0,3)) = 1;
	grid.isonb(grid.nv(grid.nbe(:,3) == 0,1)) = 1;
	grid.isonb(grid.nv(grid.nbe(:,3) == 0,2)) = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Find the closest ## elements to each element and save them.              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stridx=find(grid.ncfile=='/');
grid.name=grid.ncfile(max(stridx)+1:length(grid.ncfile)-8);

if (exist(['grid_closest/' grid.name '.mat'],'file')==0)
    disp('Calculating closest 40 elements to each element. This will take time for large n. Scale as O(n^2*log(n)) roughly.')

    nbr_closest_elem = 40;
    grid.closest=zeros(grid.nele,nbr_closest_elem);	
tic
    for i=1:grid.nele
			dist=(grid.xc-grid.xc(i)).^2 + (grid.yc-grid.yc(i)).^2;
            [~, idx]=sort(dist,'ascend');
            grid.closest(i,:) = idx(2:(nbr_closest_elem+1));
	end
toc
    
    gridclosest=grid.closest;
    if ~exist('grid_closest','dir'), mkdir('grid_closest'); end
    save(['grid_closest/' grid.name '.mat'],'gridclosest');

else
	disp('Loading existing file with closest element data')
    load(['grid_closest/' grid.name '.mat']);
    grid.closest=gridclosest;

end


