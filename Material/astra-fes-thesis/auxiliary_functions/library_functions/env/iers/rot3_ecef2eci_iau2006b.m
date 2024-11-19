function ROT = rot3_ecef2eci_iau2006b(tt_s, ut1_s, xp, yp)

    % Compute the rotation matrix from the `ITRF` (i.e., ECEF) to the 
    % `GCRF` (i.e., ECI) according to the IAU2006/2000B model.
    % The IAU2006B model requires the interpolated values of the pole coordinates
    % xp and yp from the latest EOP data. It neglects the Free Core Nutation (FCN)
    % corrections (dX and dY). The precession-nutation matrix is computed following 
    % the IAU 2006A model but with truncated expressions for the nutation corrections. 
    % It is slighlty more precise than the original IAU2000B model across the 1995-2050 
    % time-span.
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
    % Notes
    % -----
    %   The nutation series for the IAU 2000B model is truncated from nearly 1400 terms 
    %   to only 77, yet it still delivers results of 1 mas accuracy at present epochs. 
    %   In particular, it delivers a pole accurate to 1 mas from 1990 to 2100 (only very 
    %   occasionally just outside 1 mas). The coefficients are taken from SOFA's implementation, 
    %   which slightly differs from those reported in McCarthy and Luzum (2003). Comparisons 
    %   with the IAU 2006A show that the SOFA version between 1995 and 2050 delivers 0.283 
    %   mas RMSE (0.994 mas in the worst case), whereas the IERS conventions website version 
    %   delivers 0.312 RMSE (1.125 mas in the worst case). 
    %
    %   A simplified version of the Fundamental arguments, taken from Simon et al (1994) is 
    %   exploited as the error introduced is below the model accuracy (~0.1 mas).
    %
    % References 
	% ----------
	%   [1] Luzum, B. and Petit G. (2012), The IERS Conventions (2010)
    %   [2] Simon, J. et al. (1994), Numerical expressions for precession formulae and mean 
    %       elements for the Moon and the planets.
    %   [3] Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent 
    %       with the IAU 2006 resolutions. 
    
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

    % Compute CIP coordinates
    [x, y] = cip_coords(t); 

    % Compute CIO Locator
    s = cio_locator(t, x, y);

    % Form intermediate-to-celestial matrix 
    Q = xys2m(x, y, s); 

end

function [x, y] = cip_coords(t)

    % Compute the CIP coordinates for the CPNd model. Time `t` is 
    % expressed in TT Julian centuries since J2000.

    % Computes Fukushima-Williams angles
    [gam, phi, psi, obl] = fw_angles(t);
    
    % Computes IAU 2000 nutation components
    [dpsi, dobl] = nutation(t);

    % Retrieves CIP coordinates by applying IAU-2006 compatible nutations 
    [x, y] = fw2xy(gam, phi, psi + dpsi, obl + dobl);

end

function [dpsi, dobl] = nutation(t)
   
    % Compute luni-solar nutation contributions
    [dpsi_ls, dobl_ls] = nutation00(t);
    
    % Add offset to account for truncated planetary contributions
    dpsi_pl = arcsec2rad(-0.135*1e-3); 
    dobl_pl = arcsec2rad( 0.388*1e-3);
    
    % Sum contributions
    dpsi = dpsi_ls + dpsi_pl; 
    dobl = dobl_ls + dobl_pl;
    
end

function s = cio_locator(t, x, y)
    
    % Compute Fundamental arguments (uses 2006A arguments)
    fa = fundamental_arguments_a(t); 
    
    % Retrieve fundamental arguments used
    Ma = fa(1); 
    Sa = fa(2);
    um = fa(3); 
    Ds = fa(4); 
    Om = fa(5); 
    
    l_Ve = fa(7); 
    l_Ea = fa(8);
    pa = fa(14);
    
    Ma2 = 2 * Ma;
    um1 = 2 * um;
    um2 = 4 * um;
    um4 = -2 * um;
    Ds1 = -2 * Ds;
    Ds2 = -4 * Ds;
    Ds3 = -1 * Ds;
    Ds4 = 2 * Ds;
    Om2 = 2 * Om;
    Om3 = 3 * Om;
    Om4 = -1 * Om;
    Om5 = 4 * Om;
    Om6 = -3 * Om;
    Om7 = -2 * Om;

    l_Ve1 = -8 * l_Ve;
    l_Ve2 = 8 * l_Ve;
    l_Ea1 = 12 * l_Ea;
    l_Ea2 = -13 * l_Ea;

    pa1 = -1 * pa;
    x1x1 = 94.0;
    x1x2 = 3808.35;
    x1x3 = -119.94;
    x1x4 = -72574.09;
    x1x5 = 27.7;
    x1x6 = 15.61;
    arg = Om;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 - 2640.73*sarg + 0.39*carg;
    x1x2 = x1x2 + 1.71*sarg - 0.03*carg;
    x1x3 = x1x3 + 743.53*sarg - 0.17*carg;
    x1x4 = x1x4 + 0.3*sarg -23.51*carg;
    x1x5 = x1x5 - 0.26*sarg - 0.01*carg;
    arg = Om2;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 - 63.53*sarg + 0.02*carg;
    x1x2 = x1x2 - 0.07*sarg + 3.57*carg;
    x1x3 = x1x3 - 8.85*sarg + 0.01*carg;
    x1x4 = x1x4 + (0.22*carg);
    arg = um1 + Ds1 + Om3;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 - 11.75*sarg - 0.01*carg;
    x1x2 = x1x2 + (0.48*carg);
    arg = um1 + Ds1 + Om;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 - 11.21*sarg - 0.01*carg;
    x1x3 = x1x3 + (- 0.55*sarg);
    arg = um1 + Ds1 + Om2;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + (4.57*sarg);
    x1x3 = x1x3 + 56.91*sarg + 0.06*carg;
    x1x4 = x1x4 - 0.03*sarg -1.39*carg;
    arg = um1 + Om3;
    sarg = sin(arg);
    x1x1 = x1x1 + (-2.02*sarg);
    arg = um1 + Om;
    sarg = sin(arg);
    x1x1 = x1x1 + (-1.98*sarg);
    x1x3 = x1x3 + (1.67*sarg);
    arg = Om3;
    sarg = sin(arg);
    x1x1 = x1x1 + (1.72*sarg);
    arg = Sa + Om;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + 1.41*sarg + 0.01*carg;
    arg = Sa + Om4;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + 1.26*sarg + 0.01*carg;
    arg = Ma + Om4;
    sarg = sin(arg);
    x1x1 = x1x1 + (0.63*sarg);
    x1x3 = x1x3 + (- 0.25*sarg);
    arg = Ma + Om;
    sarg = sin(arg);
    x1x1 = x1x1 + (0.63*sarg);
    x1x3 = x1x3 + (- 0.27*sarg);
    arg = Sa + um1 + Ds1 + Om3;
    sarg = sin(arg);
    x1x1 = x1x1 + (- 0.46*sarg);
    arg = Sa + um1 + Ds1 + Om;
    sarg = sin(arg);
    x1x1 = x1x1 + (- 0.45*sarg);
    arg = um2 + Ds2 + Om5;
    sarg = sin(arg);
    x1x1 = x1x1 + (- 0.36*sarg);
    arg = um + Ds3 + Om + l_Ve1 + l_Ea1;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + 0.24*sarg + 0.12*carg;
    arg = um1;
    sarg = sin(arg);
    x1x1 = x1x1 + (- 0.32*sarg);
    x1x3 = x1x3 + (- 0.11*sarg);
    arg = um1 + Om2;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + (- 0.28*sarg);
    x1x3 = x1x3 + 9.84*sarg - 0.01*carg;
    x1x4 = x1x4 - 0.01*sarg - 0.24*carg;
    arg = Ma + um1 + Om3;
    sarg = sin(arg);
    x1x1 = x1x1 + (- 0.27*sarg);
    arg = Ma + um1 + Om;
    sarg = sin(arg);
    x1x1 = x1x1 + (- 0.26*sarg);
    x1x3 = x1x3 + (0.22*sarg);
    arg = um1 + Ds1;
    sarg = sin(arg);
    x1x1 = x1x1 + (0.21*sarg);
    arg = Sa + um4 + Ds4 + Om6;
    sarg = sin(arg);
    x1x1 = x1x1 + (- 0.19*sarg);
    arg = Sa + um4 + Ds4 + Om4;
    sarg = sin(arg);
    x1x1 = x1x1 + (- 0.18*sarg);
    arg = l_Ve2 + l_Ea2 + pa1;
    carg = cos(arg);
    sarg = sin(arg);
    x1x1 = x1x1 + 0.1*sarg - 0.05*carg;
    arg = Ds4;
    sarg = sin(arg);
    x1x1 = x1x1 + (- 0.15*sarg);
    x1x3 = x1x3 + (- 0.27*sarg);
    arg = Ma2 + um4 + Om4;
    sarg = sin(arg);
    x1x1 = x1x1 + (0.14*sarg);
    x1x3 = x1x3 + (0.2*sarg);
    arg = Sa + um1 + Ds1 + Om2;
    sarg = sin(arg);
    x1x1 = x1x1 + (0.14*sarg);
    x1x3 = x1x3 + (2.23*sarg);
    arg = Ma + Ds1 + Om;
    sarg = sin(arg);
    x1x1 = x1x1 + (- 0.14*sarg);
    arg = Ma + Ds1 + Om4;
    sarg = sin(arg);
    x1x1 = x1x1 + (- 0.14*sarg);
    arg = um2 + Ds1 + Om5;
    sarg = sin(arg);
    x1x1 = x1x1 + (- 0.13*sarg);
    arg = um1 + Ds1 + Om5;
    sarg = sin(arg);
    x1x1 = x1x1 + (0.11*sarg);
    arg = Ma + um4 + Om6;
    sarg = sin(arg);
    x1x1 = x1x1 + (- 0.11*sarg);
    arg = Ma + um4 + Om4;
    sarg = sin(arg);
    x1x1 = x1x1 + (- 0.11*sarg);
    arg = Sa;
    carg = cos(arg);
    sarg = sin(arg);
    x1x3 = x1x3 - 6.38*sarg - 0.05*carg;
    arg = Ma;
    sarg = sin(arg);
    x1x3 = x1x3 + (-3.07*sarg);
    arg = Ma + um1 + Om2;
    sarg = sin(arg);
    x1x3 = x1x3 + (1.3*sarg);
    arg = Sa + um4 + Ds4 + Om7;
    sarg = sin(arg);
    x1x3 = x1x3 + (0.93*sarg);
    arg = Ma + Ds1;
    sarg = sin(arg);
    x1x3 = x1x3 + (0.68*sarg);
    arg = Ma + um4 + Om7;
    sarg = sin(arg);
    x1x3 = x1x3 + (0.53*sarg);
    arg = Ma + um4 + Ds1 + Om7;
    sarg = sin(arg);
    x1x3 = x1x3 + (- 0.26*sarg);
    arg = Ma2 + Ds1;
    sarg = sin(arg);
    x1x3 = x1x3 + (- 0.21*sarg);
    arg = um1 + Ds4 + Om2;
    sarg = sin(arg);
    x1x3 = x1x3 + (0.17*sarg);
    arg = Ma2 + um1 + Om2;
    sarg = sin(arg);
    x1x3 = x1x3 + (0.13*sarg);
    arg = Ma2;
    sarg = sin(arg);
    x1x3 = x1x3 + (- 0.13*sarg);
    arg = Ma + um1 + Ds1 + Om2;
    sarg = sin(arg);
    x1x3 = x1x3 + (- 0.12*sarg);

    s = x1x1 + t*(x1x2 + t*(x1x3 + t*(x1x4 + t*(x1x5 + t*x1x6))));

    % Micro arcseconds to radians conversion factor
    muas2rad = 1e-6*pi/648000;

    s = s * muas2rad - 0.5*x*y;
end


function [dpsi, dobl] = nutation00(t)

    [Ma, Sa, um, Ds, Om] = fundamental_arguments_b(t);
    
    Ma2 = -1 * Ma;
    Ma3 = -2 * Ma;
    Ma4 = 2 * Ma;
    Ma5 = 3 * Ma;
    Sa2 = -1 * Sa;
    Sa3 = -2 * Sa;
    Sa4 = 2 * Sa;
    um1 = 2 * um;
    um2 = -2 * um;
    Ds1 = -2 * Ds;
    Ds2 = 2 * Ds;
    Ds4 = 4 * Ds;
    Om2 = 2 * Om;

    dpsi0 = 0.0;
    dobl0 = 0.0;
    dpsi1 = 0.0;
    dobl1 = 0.0;
    
    arg = Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 1.72064161e7*sarg + 3338.6*carg;
    dobl0 = dobl0 + 1537.7*sarg + 9.2052331e6*carg;
    dpsi1 = dpsi1 + (-17466.6*sarg);
    dobl1 = dobl1 + (908.6*carg);
    arg = um1 + Ds1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 1.3170906e6*sarg - 1369.6*carg;
    dobl0 = dobl0 - 458.7*sarg + 573033.6*carg;
    dpsi1 = dpsi1 + (-167.5*sarg);
    dobl1 = dobl1 + (-301.5*carg);
    arg = um1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 227641.3*sarg + 279.6*carg;
    dobl0 = dobl0 + 137.4*sarg + 97845.9*carg;
    dpsi1 = dpsi1 + (-23.4*sarg);
    dobl1 = dobl1 + (-48.5*carg);
    arg = Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 207455.4*sarg - 69.8*carg;
    dobl0 = dobl0 - 29.1*sarg - 89749.2*carg;
    dpsi1 = dpsi1 + (20.7*sarg);
    dobl1 = dobl1 + (47.0*carg);
    arg = Sa;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 147587.7*sarg + 1181.7*carg;
    dobl0 = dobl0 - 192.4*sarg + 7387.1*carg;
    dpsi1 = dpsi1 + (-363.3*sarg);
    dobl1 = dobl1 + (-18.4*carg);
    arg = Sa + um1 + Ds1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 51682.1*sarg - 52.4*carg;
    dobl0 = dobl0 - 17.4*sarg + 22438.6*carg;
    dpsi1 = dpsi1 + (122.6*sarg);
    dobl1 = dobl1 + (-67.7*carg);
    arg = Ma;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 71115.9*sarg - 87.2*carg;
    dobl0 = dobl0 + 35.8*sarg - 675.0*carg;
    dpsi1 = dpsi1 + (7.3*sarg);
    arg = um1 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 38729.8*sarg + 38.0*carg;
    dobl0 = dobl0 + 31.8*sarg + 20072.8*carg;
    dpsi1 = dpsi1 + (-36.7*sarg);
    dobl1 = dobl1 + (1.8*carg);
    arg = Ma + um1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 30146.1*sarg + 81.6*carg;
    dobl0 = dobl0 + 36.7*sarg + 12902.5*carg;
    dpsi1 = dpsi1 + (-3.6*sarg);
    dobl1 = dobl1 + (-6.3*carg);
    arg = Sa2 + um1 + Ds1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 21582.9*sarg + 11.1*carg;
    dobl0 = dobl0 + 13.2*sarg - 9592.9*carg;
    dpsi1 = dpsi1 + (-49.4*sarg);
    dobl1 = dobl1 + (29.9*carg);
    arg = um1 + Ds1 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 12822.7*sarg + 18.1*carg;
    dobl0 = dobl0 + 3.9*sarg - 6898.2*carg;
    dpsi1 = dpsi1 + (13.7*sarg);
    dobl1 = dobl1 + (-0.9*carg);
    arg = Ma2 + um1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 12345.7*sarg + 1.9*carg;
    dobl0 = dobl0 - 0.4*sarg - 5331.1*carg;
    dpsi1 = dpsi1 + (1.1*sarg);
    dobl1 = dobl1 + (3.2*carg);
    arg = Ma2 + Ds2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 15699.4*sarg - 16.8*carg;
    dobl0 = dobl0 + 8.2*sarg - 123.5*carg;
    dpsi1 = dpsi1 + (1.0*sarg);
    arg = Ma + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 6311.0*sarg + 2.7*carg;
    dobl0 = dobl0 - 0.9*sarg - 3322.8*carg;
    dpsi1 = dpsi1 + (6.3*sarg);
    arg = Ma2 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 5797.6*sarg - 18.9*carg;
    dobl0 = dobl0 - 7.5*sarg + 3142.9*carg;
    dpsi1 = dpsi1 + (-6.3*sarg);
    arg = Ma2 + um1 + Ds2 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 5964.1*sarg + 14.9*carg;
    dobl0 = dobl0 + 6.6*sarg + 2554.3*carg;
    dpsi1 = dpsi1 + (-1.1*sarg);
    dobl1 = dobl1 + (-1.1*carg);
    arg = Ma + um1 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 5161.3*sarg + 12.9*carg;
    dobl0 = dobl0 + 7.8*sarg + 2636.6*carg;
    dpsi1 = dpsi1 + (-4.2*sarg);
    arg = Ma3 + um1 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 4589.3*sarg + 3.1*carg;
    dobl0 = dobl0 + 2.0*sarg - 2423.6*carg;
    dpsi1 = dpsi1 + (5.0*sarg);
    dobl1 = dobl1 + (-1.0*carg);
    arg = Ds2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 6338.4*sarg - 15.0*carg;
    dobl0 = dobl0 + 2.9*sarg - 122.0*carg;
    dpsi1 = dpsi1 + (1.1*sarg);
    arg = um1 + Ds2 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 3857.1*sarg + 15.8*carg;
    dobl0 = dobl0 + 6.8*sarg + 1645.2*carg;
    dpsi1 = dpsi1 + (-0.1*sarg);
    dobl1 = dobl1 + (-1.1*carg);
    arg = Sa3 + um1 + Ds1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + (3248.1*sarg);
    dobl0 = dobl0 + (-1387.0*carg);
    arg = Ma3 + Ds2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 4772.2*sarg - 1.8*carg;
    dobl0 = dobl0 - 2.5*sarg + 47.7*carg;
    arg = Ma4 + um1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 3104.6*sarg + 13.1*carg;
    dobl0 = dobl0 + 5.9*sarg + 1323.8*carg;
    dpsi1 = dpsi1 + (-0.1*sarg);
    dobl1 = dobl1 + (-1.1*carg);
    arg = Ma + um1 + Ds1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 2859.3*sarg - 0.1*carg;
    dobl0 = dobl0 - 0.3*sarg - 1233.8*carg;
    dobl1 = dobl1 + (1.0*carg);
    arg = Ma2 + um1 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 2044.1*sarg + 1.0*carg;
    dobl0 = dobl0 - 0.3*sarg - 1075.8*carg;
    dpsi1 = dpsi1 + (2.1*sarg);
    arg = Ma4;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 2924.3*sarg - 7.4*carg;
    dobl0 = dobl0 + 1.3*sarg - 60.9*carg;
    arg = um1;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 2588.7*sarg - 6.6*carg;
    dobl0 = dobl0 + 1.1*sarg - 55.0*carg;
    arg = Sa + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 1405.3*sarg + 7.9*carg;
    dobl0 = dobl0 - 4.5*sarg + 855.1*carg;
    dpsi1 = dpsi1 + (-2.5*sarg);
    dobl1 = dobl1 + (-0.2*carg);
    arg = Ma2 + Ds2 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 1516.4*sarg + 1.1*carg;
    dobl0 = dobl0 - 0.1*sarg - 800.1*carg;
    dpsi1 = dpsi1 + (1.0*sarg);
    arg = Sa4 + um1 + Ds1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 1579.4*sarg - 1.6*carg;
    dobl0 = dobl0 - 0.5*sarg + 685.0*carg;
    dpsi1 = dpsi1 + (7.2*sarg);
    dobl1 = dobl1 + (-4.2*carg);
    arg = um2 + Ds2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 2178.3*sarg + 1.3*carg;
    dobl0 = dobl0 + 1.3*sarg - 16.7*carg;
    arg = Ma + Ds1 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 1287.3*sarg - 3.7*carg;
    dobl0 = dobl0 - 1.4*sarg + 695.3*carg;
    dpsi1 = dpsi1 + (-1.0*sarg);
    arg = Sa2 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 1265.4*sarg + 6.3*carg;
    dobl0 = dobl0 + 2.6*sarg + 641.5*carg;
    dpsi1 = dpsi1 + (1.1*sarg);
    arg = Ma2 + um1 + Ds2 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 1020.4*sarg + 2.5*carg;
    dobl0 = dobl0 + 1.5*sarg + 522.2*carg;
    arg = Sa4;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 1670.7*sarg - 1.0*carg;
    dobl0 = dobl0 + 1.0*sarg + 16.8*carg;
    dpsi1 = dpsi1 + (-8.5*sarg);
    dobl1 = dobl1 + (-0.1*carg);
    arg = Ma + um1 + Ds2 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 769.1*sarg + 4.4*carg;
    dobl0 = dobl0 + 1.9*sarg + 326.8*carg;
    arg = Ma3 + um1;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 1102.4*sarg - 1.4*carg;
    dobl0 = dobl0 + 0.2*sarg + 10.4*carg;
    arg = Sa + um1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 756.6*sarg - 1.1*carg;
    dobl0 = dobl0 - 0.5*sarg - 325.0*carg;
    dpsi1 = dpsi1 + (-2.1*sarg);
    arg = um1 + Ds2 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 663.7*sarg + 2.5*carg;
    dobl0 = dobl0 + 1.4*sarg + 335.3*carg;
    dpsi1 = dpsi1 + (-1.1*sarg);
    arg = Sa2 + um1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 714.1*sarg + 0.8*carg;
    dobl0 = dobl0 + 0.4*sarg + 307.0*carg;
    dpsi1 = dpsi1 + (2.1*sarg);
    arg = Ds2 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 630.2*sarg + 0.2*carg;
    dobl0 = dobl0 + 0.4*sarg + 327.2*carg;
    dpsi1 = dpsi1 + (-1.1*sarg);
    arg = Ma + um1 + Ds1 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 580.0*sarg + 0.2*carg;
    dobl0 = dobl0 - 0.1*sarg - 304.5*carg;
    dpsi1 = dpsi1 + (1.0*sarg);
    arg = Ma4 + um1 + Ds1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 644.3*sarg - 0.7*carg;
    dobl0 = dobl0 - 0.4*sarg - 276.8*carg;
    arg = Ma3 + Ds2 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 577.4*sarg - 1.5*carg;
    dobl0 = dobl0 - 0.5*sarg + 304.1*carg;
    dpsi1 = dpsi1 + (-1.1*sarg);
    arg = Ma4 + um1 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 535.0*sarg + 2.1*carg;
    dobl0 = dobl0 + 1.2*sarg + 269.5*carg;
    arg = Sa2 + um1 + Ds1 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 475.2*sarg - 0.3*carg;
    dobl0 = dobl0 - 0.3*sarg + 271.9*carg;
    dpsi1 = dpsi1 + (-1.1*sarg);
    arg = Ds1 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 494.0*sarg - 2.1*carg;
    dobl0 = dobl0 - 0.9*sarg + 272.0*carg;
    dpsi1 = dpsi1 + (-1.1*sarg);
    arg = Ma2 + Sa2 + Ds2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 735.0*sarg - 0.8*carg;
    dobl0 = dobl0 + 0.4*sarg - 5.1*carg;
    arg = Ma4 + Ds1 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 406.5*sarg + 0.6*carg;
    dobl0 = dobl0 + 0.1*sarg - 220.6*carg;
    arg = Ma + Ds2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 657.9*sarg - 2.4*carg;
    dobl0 = dobl0 + 0.2*sarg - 19.9*carg;
    arg = Sa + um1 + Ds1 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 357.9*sarg + 0.5*carg;
    dobl0 = dobl0 + 0.1*sarg - 190.0*carg;
    arg = Ma + Sa2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 472.5*sarg - 0.6*carg;
    dobl0 = dobl0 + 0.3*sarg - 4.1*carg;
    arg = Ma3 + um1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 307.5*sarg - 0.2*carg;
    dobl0 = dobl0 - 0.1*sarg + 131.3*carg;
    arg = Ma5 + um1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 290.4*sarg + 1.5*carg;
    dobl0 = dobl0 + 0.7*sarg + 123.3*carg;
    arg = Sa2 + Ds2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 434.8*sarg - 1.0*carg;
    dobl0 = dobl0 + 0.2*sarg - 8.1*carg;
    arg = Ma + Sa2 + um1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 287.8*sarg + 0.8*carg;
    dobl0 = dobl0 + 0.4*sarg + 123.2*carg;
    arg = Ds;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 423.0*sarg + 0.5*carg;
    dobl0 = dobl0 - 0.2*sarg - 2.0*carg;
    arg = Ma2 + Sa2 + um1 + Ds2 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 281.9*sarg + 0.7*carg;
    dobl0 = dobl0 + 0.3*sarg + 120.7*carg;
    arg = Ma2 + um1;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 405.6*sarg + 0.5*carg;
    dobl0 = dobl0 - 0.2*sarg + 4.0*carg;
    arg = Sa2 + um1 + Ds2 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 264.7*sarg + 1.1*carg;
    dobl0 = dobl0 + 0.5*sarg + 112.9*carg;
    arg = Ma3 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 229.4*sarg - 1.0*carg;
    dobl0 = dobl0 - 0.4*sarg + 126.6*carg;
    arg = Ma + Sa + um1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 248.1*sarg - 0.7*carg;
    dobl0 = dobl0 - 0.3*sarg - 106.2*carg;
    arg = Ma4 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 217.9*sarg - 0.2*carg;
    dobl0 = dobl0 - 0.2*sarg - 112.9*carg;
    arg = Ma2 + Sa + Ds;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 327.6*sarg + 0.1*carg;
    dobl0 = dobl0 + (-0.9*carg);
    arg = Ma + Sa;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 338.9*sarg + 0.5*carg;
    dobl0 = dobl0 - 0.2*sarg + 3.5*carg;
    arg = Ma + um1;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 333.9*sarg - 1.3*carg;
    dobl0 = dobl0 + 0.1*sarg - 10.7*carg;
    arg = Ma2 + um1 + Ds1 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 198.7*sarg - 0.6*carg;
    dobl0 = dobl0 - 0.2*sarg + 107.3*carg;
    arg = Ma + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + (-198.1*sarg);
    dobl0 = dobl0 + (85.4*carg);
    arg = Ma2 + Ds;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 402.6*sarg - 35.3*carg;
    dobl0 = dobl0 - 13.9*sarg - 55.3*carg;
    arg = um1 + Ds + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 166.0*sarg - 0.5*carg;
    dobl0 = dobl0 - 0.2*sarg - 71.0*carg;
    arg = Ma2 + um1 + Ds4 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 152.1*sarg + 0.9*carg;
    dobl0 = dobl0 + 0.4*sarg + 64.7*carg;
    arg = Ma2 + Sa + Ds + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + (131.4*sarg);
    dobl0 = dobl0 + (-70.0*carg);
    arg = Sa3 + um1 + Ds1 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + (-128.3*sarg);
    dobl0 = dobl0 + (67.2*carg);
    arg = Ma + um1 + Ds2 + Om;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 - 133.1*sarg + 0.8*carg;
    dobl0 = dobl0 + 0.4*sarg + 66.3*carg;
    arg = Ma3 + um1 + Ds2 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 138.3*sarg - 0.2*carg;
    dobl0 = dobl0 - 0.2*sarg - 59.4*carg;
    arg = Ma2 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + 140.5*sarg + 0.4*carg;
    dobl0 = dobl0 + 0.2*sarg - 61.0*carg;
    arg = Ma + Sa + um1 + Ds1 + Om2;
    carg = cos(arg); 
    sarg = sin(arg); 
    dpsi0 = dpsi0 + (129.0*sarg);
    dobl0 = dobl0 + (-55.6*carg);

    % Conversion factor
    muas2rad = 1e-6*pi/648000;     

    dpsi = (dpsi0 + t*dpsi1) * muas2rad; 
    dobl = (dobl0 + t*dobl1) * muas2rad; 

end