function W = pm_rotm(t_c, xp, yp)

    % Compute the Polar Motion rotation matrix from ITRF to TIRS according
    % to the IAU 2010 Conventions. 
    
    % Parameters 
    % ----------
    %   t_c: double
    %       Time expressed as TT Julian centuries since J2000
    %   xp: double
    %       Polar x coordinate [rad]
    %   yp: double 
    %       Polar y coordinate [rad]
    % 
    % Returns 
    % -------
    %   W: double (3, 3) 
    %       Polar motion rotation matrix
    %
    % References 
	% ----------
	% [1] Luzum, B. and Petit G. (2012), The IERS Conventions (2010)
    
    sp = arcsec2rad(-47e-6*t_c);
    W = angle3_to_dcm(yp, xp, -sp, "XYZ");
    
end