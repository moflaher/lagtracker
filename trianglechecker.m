function [isintriangle]=trianglechecker(xt,yt,x0,y0) 
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: xt,yt,x0,y0
% Return: triangle
%
%   Check to see if x0,y0 is in element define by xt,yt.
%   Using Algorithm Used for Scene Rendering in Computer Graphics   
%
%   This should be speed tested for current use cases against buildin matlab functions.
%   Very similar to other triangle checkers should be a way to merge them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%code taken from pt_farm changed syntax for matlab
	%==============================================================================|
  	%  Determine if Point (x0,y0) Is In triangle defined by nodes (xt(3),yt(3))    |
  	%  Using Algorithm Used for Scene Rendering in Computer Graphics               |
  	%==============================================================================|
  	%------------------------------------------------------------------------------|


  	isintriangle = 0; 

  	f1 = (y0-yt(1))*(xt(2)-xt(1)) - (x0-xt(1))*(yt(2)-yt(1));
  	f2 = (y0-yt(3))*(xt(1)-xt(3)) - (x0-xt(3))*(yt(1)-yt(3));
  	f3 = (y0-yt(2))*(xt(3)-xt(2)) - (x0-xt(2))*(yt(3)-yt(2));

  	if(f1*f3 >= 0 && f3*f2 >= 0)
		isintriangle = 1;
	end



