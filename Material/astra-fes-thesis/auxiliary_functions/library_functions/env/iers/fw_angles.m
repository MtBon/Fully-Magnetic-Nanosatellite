function [gam, phi, psi, obl] = fw_angles(t)

	% Compute the precession angles in radians, following the IAU 2006 
	% Fukushima-Williams formulation. 
	%
	% Parameters
	% ----------
	% 	t: double
	%		Time expressed in TT centuries since J2000
	%
	% Returns 
	% -------
	% 	gam: double
	%		FW 1st angle [rad]
	%	phi: double
	%		FW 2nd angle [rad]
	%	psi: double
	%		FW 3rd angle [rad]
	%	obl: double
	%		FW 4th angle [rad]
	%
	% References 
	% ----------
	% [1] Luzum, B. and Petit G. (2012), The IERS Conventions (2010)
	% [2] Wallace P. T. and Capitaine N. (2006), Precession-nutation
	% 	procedures with IAU 20006 resolutions
	
    gam = -0.052928 + t*(10.556378 + t*(0.4932044 + t*(-0.00031238 + ...
          t*(-0.000002788 + 0.0000000260*t))));

    phi = 84381.412819 + t*(-46.811016 + t*(0.0511268 + t*(0.00053289 - ...
          t*(0.000000440 + 0.0000000176*t))));

    psi = -0.041775 + t*(5038.481484 + t*(1.5584175 - t*(0.00018522 + ...
          t*(0.000026452 + 0.0000000148*t))));
    
    obl = 84381.406 + t*(-46.836769 + t*(-0.0001831 + t*(0.00200340 - ...
          t*(0.000000576 + 0.0000000434*t))));

    % Conver to radians
    arcsec2rad = pi/648000;
    
    gam = gam*arcsec2rad;
    phi = phi*arcsec2rad;
    psi = psi*arcsec2rad;
    obl = obl*arcsec2rad;
      
end