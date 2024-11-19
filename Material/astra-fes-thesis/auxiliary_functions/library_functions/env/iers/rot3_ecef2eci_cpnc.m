function ROT = rot3_ecef2eci_cpnc(tt_s, ut1_s, xp, yp)

    % Compute the rotation matrix from the `ITRF` (i.e., ECEF) to the 
    % `GCRF` (i.e., ECI) according to the CPNc approximate model.
    % The CPNc version is a concise model with a cut-off at 2.5 mas of the X 
    % and Y series, delivering a worst-case accuracy of about 15 mas between 1995-2050. 
    % It does not take into account the Free Core Nutation (~0.2 mas)
    %
    % Parameters
    % ----------
    %   tt_s: double
    %       Time expressed in TT seconds since J2000
    %   ut1_s: double
    %       Time expressed in UT1 seconds since J2000
    %   xp: double
    %       IERS's EOP interpolated Pole x-coordinate [rad]
    %   yp: double
    %       IERS's EOP interpolated Pole y-coordinate [rad]
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
    
    W = pm_rotm(tt_c, xp, yp);
    R = era_rotm(ut1_d); 
    Q = cip_motion(tt_c); 
    
    ROT = Q*R*W;
   
end

function Q = cip_motion(t)

    % Compute the CIRS to GCRS rotation matrix, according to the CPNd model for the
    % CIP coordinates. The time is expressed in TT Julian centuries since J2000.

    [x, y] = cip_coords(t); 
    s = cio_locator(t, x, y);
    Q = xys2m(x, y, s); 

end


function s = cio_locator(t, x, y)

    % Compute the CIO locator `s` in radians, given the CIP coordinates
    % `X` and `Y` and time `t' expressed in TT Julian centuries since J2000

    W =  2.182439196616 - 33.7570459536*t;

    s = t*(3809 - 72574.0*t^2) - 2641*sin(W);

    % Transform s from micro arcsec to radians 
    s = s*1e-6*pi/648000; 

    % Subtract XY/2 term from the series
    s = s - 0.5*x*y;

end

function [x, y] = cip_coords(t)

    % Compute the CIP coordinates for the CPNc model in radians. 
    % Time `t` is expressed in TT Julian centuries since J2000.

    [Ma, Sa, um, Ds, Wm] = fundamental_arguments_b(t);

    um1 = 2 * um;
    um2 = -2 * um;
    Ds1 = 2 * Ds;
    Ds2 = -2 * Ds;
    Wm2 = 2 * Wm;
    Wm3 = -2 * Wm;
    Wm4 = -1 * Wm;

    x1x1 = -17251.0;
    x1x2 = 2.004191898e9;
    x1x3 = -429783.0;
    x1x4 = -198618.0;
    x2x1 = -5530.0;
    x2x2 = -25896.0;
    x2x3 = -2.2407275e7;
    x2x4 = 1901.0;
    x2x5 = 1113.0;
    
    arg = Wm;
    carg = cos(arg); 
    sarg = sin(arg);
    x1x1 = x1x1 + (-6.844318e6*sarg);
    x1x2 = x1x2 + -3310.0*sarg + 205833.0*carg;
    x2x1 = x2x1 + (9.205236e6*carg);
    x2x2 = x2x2 + (153042.0*sarg);
    arg = Wm2;
    carg = cos(arg); 
    sarg = sin(arg);
    x1x1 = x1x1 + (82169.0*sarg);
    x2x1 = x2x1 + (-89618.0*carg);
    arg = Ds1;
    sarg = sin(arg);
    x1x1 = x1x1 + (2521.0*sarg);
    arg = um1 + Ds2 + Wm;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + (5096.0*sarg);
    x2x1 = x2x1 + (-6918.0*carg);
    arg = um1 + Ds2 + Wm2;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + (-523908.0*sarg);
    x1x2 = x1x2 + (12814.0*carg);
    x2x1 = x2x1 + (573033.0*carg);
    x2x2 = x2x2 + (11714.0*sarg);
    arg = um1 + Wm;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + (-15407.0*sarg);
    x2x1 = x2x1 + (20070.0*carg);
    arg = um1 + Wm2;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + (-90552.0*sarg);
    x2x1 = x2x1 + (97847.0*carg);
    arg = Sa + um2 + Ds1 + Wm3;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + (-8585.0*sarg);
    x2x1 = x2x1 + (-9593.0*carg);
    arg = Sa;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + (58707.0*sarg);
    x2x1 = x2x1 + (7387.0*carg);
    arg = Sa + um1 + Ds2 + Wm2;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + (-20558.0*sarg);
    x2x1 = x2x1 + (22438.0*carg);
    arg = Ma + um2 + Wm3;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + (-4911.0*sarg);
    x2x1 = x2x1 + (-5331.0*carg);
    arg = Ma + Ds2;
    sarg = sin(arg);
    x1x1 = x1x1 + (-6245.0*sarg);
    arg = Ma;
    sarg = sin(arg);
    x1x1 = x1x1 + (28288.0*sarg);
    arg = Ma + Wm;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + (2512.0*sarg);
    x2x1 = x2x1 + (-3324.0*carg);
    arg = Ma + um1 + Wm2;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + (-11992.0*sarg);
    x2x1 = x2x1 + (12903.0*carg);
    arg = Ma + um2 + Ds2 + Wm3;
    carg = cos(arg);
    x2x1 = x2x1 + (2555.0*carg);
    arg = Ma + Wm4;
    carg = cos(arg);
    x2x1 = x2x1 + (3144.0*carg);
    arg = Ma + um1 + Wm;
    carg = cos(arg);
    x2x1 = x2x1 + (2636.0*carg);

    x = x1x1 + t*(x1x2 + t*(x1x3 + t*x1x4)); 
    y = x2x1 + t*(x2x2 + t*(x2x3+ t*(x2x4 + t*x2x5))); 
    
    % Micro arcseconds to radians conversion factor
    muas2rad = 1e-6*pi/648000;
        
    x = x * muas2rad;
    y = y * muas2rad;
    
end