function R = era_rotm(t)
    
    % Compute the TIRS to CIRS Earth Rotation matrix, according to the IERS
    % 2010 conventions at time `t'. 
    %
    % Parameters
    % ----------
    %   t: double
    %       Time expressed as UT1 days since J2000.0
    %
    % Returns 
    % -------
    %   R: double (3, 3)
    %       TIRS to CIRS Earth Rotation Matrix.
    %
    % References 
	% ----------
	% [1] Luzum, B. and Petit G. (2012), The IERS Conventions (2010)

    R = angle1_to_dcm(-era_rot_angle(t), "Z");

end