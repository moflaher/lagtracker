function [lag]=calculate_swimspeed(lag)

%!==============================================================================|
%!  Calcule the swim speed in sigma coordinates, omega_swim, from a given 
%|  swim speed in z coordinates using that
%!    omega = w %%%%%%%%%- sig*del/dt
%!==============================================================================|
  
%    lag.omega_swimspeed =lag.w_swimspeed-(1+lag.sigpt).*lag.dedtp;
    lag.omega_swimspeed =lag.w_swimspeed;

   
end


