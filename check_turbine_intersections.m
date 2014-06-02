function [lag]=check_turbine_intersections(lag,grid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: lag,grid
% Return: lag
% 
%  Finds all particles in a turbine element.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:grid.nturbines
	%check if point is in turbine triangle
	intri=zeros(grid.nele,1);
	intri=lag.host==grid.turbines(i);
	
	%check if sigma level is below top of turnbine
	lessthansigma=zeros(grid.nele,1);
	lessthansigma=lag.sigpt<lag.turbine_sigma(i);
	lag.turbine_intersects(:,i)=lag.turbine_intersects(:,i)+intri.*lessthansigma;
end
    




end



