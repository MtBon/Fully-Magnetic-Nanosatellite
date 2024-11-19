function [ran, dec, pos] = ephem_sun(time)

    % Compute the Sun position in the inertial geocentric equatorial 
    % true-of-date frame. 
    %
    % Parameters
    % ----------
    %  julian_date: double
    %       Terrestrial Time (TT) seconds since J2000.0
    % 
    % Returns
    % -------
    %  ran: double 
    %       Sun's right ascension (rad) bounded between [0, 2π]
    %  dec: double 
    %       Sun's declination (rad) bounded between [-π/2, π/2]
    %  pos: double(3, 1)
    %       Sun's position vector (km) expressed in the ECI frame
    % 
    % References
    % ----------
    % [1] T.C. Van Flandern and K.F. Pulkkinen, Low precision formulae for planetary 
    %     positions, 1979, Astrophysical Journal Supplement Series

    arcsec2rad = pi/648000;

%     DJ2000 = 2451545; % Julian Date of the Reference Epoch (J2000.0)
    djd = time/86400; % Days since J2000
    
    % Julian Centuries since 1900.0
    t = (djd / 36525) + 1; 
    
    % -- Fundamental Arguments
    % Only the fractional part of the revolution is retained to avoid the necessity 
    % of carrying additional digits into the sine and cosine functions. 
    gs = rev2rad(0.993126 + 0.00273777850*djd);
    lm = rev2rad(0.606434 + 0.03660110129*djd);
    ls = rev2rad(0.779072 + 0.00273790931*djd);
    g2 = rev2rad(0.140023 + 0.00445036173*djd);
    g4 = rev2rad(0.053856 + 0.00145561327*djd);
    g5 = rev2rad(0.056531 + 0.00023080893*djd);
    rm = rev2rad(0.347343 - 0.00014709391*djd);
    
    % Geocentric Ecliptic Longitude of the Sun (rad)
    plon = 6910 * sin(gs) + 72 * sin(2 * gs) - 17 * t * sin(gs);
    plon = plon - 7 * cos(gs - g5) + 6 * sin(lm - ls) ... 
           + 5 * sin(4 * gs - 8 * g4 + 3 * g5);
    plon = plon - 5 * cos(2 * (gs - g2)) - 4 * (sin(gs - g2) ... 
           - cos(4 * gs - 8 * g4 + 3 * g5));
    plon = plon + 3 * (sin(2 * (gs - g2)) - sin(g5) - sin(2 * (gs - g5)));

    % Add correction for nutation in longitude
    plon = ls + arcsec2rad * (plon - 17 * sin(rm));
    
    % Obliquity of the Ecliptic (rad)
    eps_0 = arcsec2rad * (84428 - 47 * t + 9 * cos(rm));

    % Geocentric distance of the Sun (km)
    rsun = 149597870.691 * (1.00014 - 0.01675 * cos(gs) - 0.00014 * cos(2 * gs));
    
    % Geocentric Equatorial Right Ascension and Declination
    a = sin(plon)*cos(eps_0);
    b = cos(plon);
       
    ran = atan2(a, b);
    dec = asin(sin(eps_0) * sin(plon));
    
    % Sun's geocentric position vector.
    pos = zeros(3, 1);

    pos(1) = rsun * cos(ran)*cos(dec);
    pos(2) = rsun * sin(ran)*cos(dec);
    pos(3) = rsun * sin(dec);
    
end