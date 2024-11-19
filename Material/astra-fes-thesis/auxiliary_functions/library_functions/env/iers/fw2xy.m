function [x, y] = fw2xy(gam, phi, psi, obl)

    % Compute the CIP X and Y coordinates from Fukushima-Williams 
    % bias-precession-nutation angles, in radians.
    %
    % Parameters
    % ----------
    %   gam: double
    %       1st FW angle [rad]
    %   phi: double
    %       2nd FW angle [rad]
    %   psi: double
    %       3rd FW angle with with IAU 2006A/B nutation corrections [rad]
    %   eps: double
    %       4th FW angle with with IAU 2006A/B nutation corrections [rad]   
    %
    % Returns 
    % -------
    %   x: double 
    %       CIP x-coordinate [rad]
    %   y: double 
    %       CIP y-coordinate [rad]
    % 
    % References 
    % ----------
	% [1] Wallace P. T. and Capitaine N. (2006), Precession-nutation
	% 	procedures with IAU 20006 resolutions
	
    sobl = sin(obl); cobl = cos(obl);
    spsi = sin(psi); cpsi = cos(psi);
    sgam = sin(gam); cgam = cos(gam);
    sphi = sin(phi); cphi = cos(phi); 

    a = (sobl*cpsi*cphi - cobl*sphi);

    x = sobl*spsi*cgam - a*sgam;
    y = sobl*spsi*sgam + a*cgam;

end