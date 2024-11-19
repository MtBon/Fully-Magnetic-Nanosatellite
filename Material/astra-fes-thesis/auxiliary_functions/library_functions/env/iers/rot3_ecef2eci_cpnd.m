function ROT = rot3_ecef2eci_cpnd(tt_s, ut1_s)

    % Compute the rotation matrix from the `ITRF` (i.e., ECEF) to the 
    % `GCRF` (i.e., ECI) according to the CPNd approximate model.
    % The CPNd model is an extremely concise formulation with an accuracy of 
    % about 1 arcsec between 1995 and 2050. It neglects polar-motion (~0.25 arcsec), 
    % the Free-Core Nutation (~0.2 mas) and the CIO locator.
    %
    % Parameters
    % ----------
    %   tt_s: double
    %       Time expressed in TT seconds since J2000
    %   ut1_s: double
    %       Time expressed in UT1 seconds since J2000 
    % 
    % Returns
    % -------
    %   ROT: double(3, 3)
    %       ECEF to ECI rotation matrix.
    %
    % References 
    % ----------
    %   [1] Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent 
    %       with the IAU 2006 resolutions. 
    %   [2] Capitaine N. and Wallace P. T. (2008), Concise CIO based precession-nutation 
    %       formulations

    % Convert TT secs since J2000 to Julian centuries 
    tt_c = tt_s/3.15576e+09;
    
    % Convert UT1 secs to days
    ut1_d = ut1_s/86400;
    
    W = pm_rotm(tt_c, 0.0, 0.0);
    R = era_rotm(ut1_d); 
    Q = cip_motion(tt_c); 
    
    ROT = Q*R*W;
   
end

function Q = cip_motion(t)

    % Compute the CIRS to GCRS rotation matrix, according to the CPNd model for the
    % CIP coordinates. The time is expressed in TT Julian centuries since J2000.

    [x, y] = cip_coords(t); 
    Q = xys2m(x, y, 0.0); 

end

function [X, Y] = cip_coords(t)

    % Compute the CIP coordinates for the CPNd model. Time `t` is 
    % expressed in TT Julian centuries since J2000.

    % Micro arcseconds to radian conversion factor
    muas2rad = 1e-6*pi/648000;
    
    % Approximated fundamental arguments as linear function of time 
    W = 2.182439196616 - 33.7570459536*t;
    A = -2.776244621014 + 1256.6639307381*t;
    
    sW = sin(W); cW = cos(W); 
    sA = sin(A); cA = cos(A); 

    % Compute approximated coordinates
    X = t*(2004191898.0 - t*(429782.9 + 198618.34*t)) - 6844318.0*sW - 523908.0*sA;
    Y = -22407275.0*t^2 + 9205236.0*cW + 573033.0*cA;

    % Apply conversion
    X = X*muas2rad; 
    Y = Y*muas2rad;

end