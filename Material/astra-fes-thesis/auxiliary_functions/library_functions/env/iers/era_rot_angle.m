function ERA = era_rot_angle(t)

    % Compute the Earth Rotation Angle (ERA), i.e., the angle between the 
    % Celeestial Intermediate Origin (CIO) and the Terrestrial Intermediate
    % Origin (TIO) at time `t`. 
    % 
    % Parameters
    % ----------
    %   t: double 
    %       Time expressed as UT1 days since J2000. 
    %
    % Returns 
    % -------
    %   ERA: double
    %       Earth Rotation Angle in [rad]. 
    % 
    % References 
    % ----------
    % Luzum B. and Petit G., (2012), The IERS Conventions (2010).

    % This function uses the fractional UT1 date to gain additional 
    % precision in the computations 
    
    f = mod(t, 1.0);
    ERA = mod(2*pi * (f + 0.7790572732640 + 0.00273781191135448*t), 2*pi);

end