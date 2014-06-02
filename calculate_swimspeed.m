function [lag]=calculate_swimspeed(lag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: lag
% Return: lag
% 
%  Calcule the swim speed in sigma coordinates, omega_swim, from a given 
%  swim speed in z coordinates using that
%    omega = w %%%%%%%%%- sig*del/dt
%
%  This is likely legacy code and will be updated shortly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%    lag.omega_swimspeed =lag.w_swimspeed-(1+lag.sigpt).*lag.dedtp;
    lag.omega_swimspeed =lag.w_swimspeed;

   
end


