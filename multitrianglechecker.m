function [triangle]=multitrianglechecker(xt,yt,x0,y0) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: xt,yt,x0,y0
% Return: triangle
%
%   Check to see which element define by xt,yt that x0,y0 is in
%   Using Algorithm Used for Scene Rendering in Computer Graphics   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    f1 = (y0-yt(:,1)).*(xt(:,2)-xt(:,1)) - (x0-xt(:,1)).*(yt(:,2)-yt(:,1));
  	f2 = (y0-yt(:,3)).*(xt(:,1)-xt(:,3)) - (x0-xt(:,3)).*(yt(:,1)-yt(:,3));
  	f3 = (y0-yt(:,2)).*(xt(:,3)-xt(:,2)) - (x0-xt(:,2)).*(yt(:,3)-yt(:,2));

  	
    triangle = ((f1.*f3 >= 0) & (f3.*f2 >= 0));
end




