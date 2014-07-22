Lagtracker is a particle tracking software.

It uses the netcdf output from FVCOM to track particle motion in a model domain.


Lagtracker is based on pt_farm a FORTRAN program designed to simulate sealice.
It was converted to Matlab in 2009 to make use of parfor to make the software parallel.

Particle motion is calculated based on advection using diffusion is optional.
Diffusion was implemented based on equations in "Lagrangian Stochastic Modeling in Coastal Oceanography" by Brickman and Smith, Jan 2002.

Comprehensive testing still needs to be done on diffusion.


