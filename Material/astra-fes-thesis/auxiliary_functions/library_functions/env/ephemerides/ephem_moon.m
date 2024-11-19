function [ran, dec, pos] = ephem_moon(time)

    % Compute the Moon position in the inertial geocentric equatorial 
    % true-of-date frame. 
    %
    % Parameters
    % ----------
    %  julian_date: double
    %         Terrestrial Time (TT) seconds since J2000.0
    % 
    % Returns
    % -------
    %  ran: double 
    %       Moon's right ascension (rad) bounded between [0, 2π]
    %  dec: double 
    %       Moon's declination (rad) bounded between [-π/2, π/2]
    %  pos: double(3, 1)
    %       Moon's position vector (km) expressed in the ECI frame
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
    gm = rev2rad(0.374897 + 0.03629164709*djd);
    fm = rev2rad(0.259091 + 0.03674819520*djd);
    em = rev2rad(0.827362 + 0.03386319198*djd);
    lm = rev2rad(0.606434 + 0.03660110129*djd);
    ls = rev2rad(0.779072 + 0.00273790931*djd);
    lv = rev2rad(0.505498 + 0.00445046867*djd);
    rm = rev2rad(0.347343 - 0.00014709391*djd);

    gm2 = 2*gm; 
    gm3 = 3*gm; 
    fm2 = 2*fm;
    em2 = 2*em;
    em4 = 2*em2;
    
    % Geocentric Ecliptic Longitude of the Moon (rad)
    l = 22640 * sin(gm) - 4586 * sin(gm - em2) + 2370 * sin(em2);
    l = l + 769 * sin(gm2) - 668 * sin(gs) - 412 * sin(fm2);
    l = l - 212 * sin(gm2 - em2) - 206 * sin(gm - em2 + gs);
    l = l + 192 * sin(gm + em2) + 165 * sin(em2 - gs);
    l = l + 148 * sin(gm - gs) - 125 * sin(em) - 110 * sin(gm + gs);
    l = l - 55 * sin(fm2 - em2) - 45 * sin(gm + fm2) + 40 * sin(gm - fm2);
    l = l - 38 * sin(gm - em4) + 36 * sin(gm3) - 31 * sin(gm2 - em4);
    l = l + 28 * sin(gm - em2 - gs) - 24 * sin(em2 + gs) + 19 * sin(gm - em);
    l = l + 18 * sin(em + gs) + 15 * sin(gm + em2 - gs) + 14 * sin(gm2 + em2);
    l = l + 14 * sin(em4) - 13 * sin(gm3 - em2) - 17 * sin(rm);
    l = l - 11 * sin(gm + 16 * ls - 18 * lv) + 10 * sin(gm2 - gs) ... 
        + 9 * sin(gm - fm2 - em2);
    l = l + 9 * (cos(gm + 16 * ls - 18 * lv) - sin(gm2 - em2 + gs)) ... 
        - 8 * sin(gm + em);
    l = l + 8 * (sin(2 * (em - gs)) - sin(gm2 + gs)) - 7 * (sin(2 * gs) ... 
        + sin(gm - 2 * (em - gs)) - sin(rm));
    l = l - 6 * (sin(gm - fm2 + em2) + sin(fm2 + em2)) ...
        - 4 * (sin(gm - em4 + gs) - t * cos(gm + 16 * ls - 18 * lv));
    l = l - 4 * (sin(gm2 + fm2) - t * sin(gm + 16 * ls - 18 * lv));
    l = l + 3 * (sin(gm - 3 * em) - sin(gm + em2 + gs) ... 
        - sin(gm2 - em4 + gs) + sin(gm - 2 * gs) + sin(gm - em2 - 2 * gs));
    l = l - 2 * (sin(gm2 - em2 - gs) + sin(fm2 - em2 + gs) - sin(gm + em4));
    l = l + 2 * (sin(4 * gm) + sin(em4 - gs) + sin(gm2 - em));
    
    plon = lm + arcsec2rad * l;

    % Geocentric Ecliptic Latitude of the Moon (rad)
    b = 18461 * sin(fm) + 1010 * sin(gm + fm) + 1000 * sin(gm - fm);
    b = b - 624 * sin(fm - em2) - 199 * sin(gm - fm - em2) ... 
        - 167 * sin(gm + fm - em2);
    b = b + 117 * sin(fm + em2) + 62 * sin(gm2 + fm) + 33 * sin(gm - fm + em2);
    b = b + 32 * sin(gm2 - fm) - 30 * sin(fm - em2 + gs) ... 
        - 16 * sin(gm2 - em2 + fm);
    b = b + 15 * sin(gm + fm + em2) + 12 * sin(fm - em2 - gs) ... 
        - 9 * sin(gm - fm - em2 + gs);
    b = b - 8 * (sin(fm + rm) - sin(fm + em2 - gs)) ... 
        - 7 * sin(gm + fm - em2 + gs);
    b = b + 7 * (sin(gm + fm - gs) - sin(gm + fm - em4));
    b = b - 6 * (sin(fm + gs) + sin(3 * fm) - sin(gm - fm - gs));
    b = b - 5 * (sin(fm + em) + sin(gm + fm + gs) + sin(gm - fm + gs) ... 
        - sin(fm - gs) - sin(fm - em));
    b = b + 4 * (sin(gm3 + fm) - sin(fm - em4)) - 3 * (sin(gm - fm - em4) ... 
        - sin(gm - 3 * fm));
    b = b - 2 * (sin(gm2 - fm - em4) + sin(3 * fm - em2) - sin(gm2 - fm + em2) ... 
        - sin(gm - fm + em2 - gs));

    plat = arcsec2rad * (b + 2 * (sin(gm2 - fm - em2) + sin(gm3 - fm)));

    % Obliquity of the Ecliptic (rad)
    eps_0 = arcsec2rad * (84428 - 47 * t + 9 * cos(rm));

    % Geocentric distance of the Moon (km)
    r = 60.36298 - 3.27746*cos(gm) - 0.57994*cos(gm - em2);
    r = r - 0.46357*cos(em2) - 0.08904*cos(gm2) + 0.03865*cos(gm2 - em2);
    r = r - 0.03237*cos(em2 - gs) - 0.02688*cos(gm + em2) - 0.02358*cos(gm - em2 + gs);
    r = r - 0.02030*cos(gm - gs) + 0.01719*cos(em) + 0.01671*cos(gm + gs);
    r = r + 0.01247*cos(gm - fm2) + 0.00704*cos(gs) + 0.00529*cos(em2 + gs);
    r = r - 0.00524*cos(gm - em4) + 0.00398*cos(gm - em2 - gs) - 0.00366*cos(gm3);
    r = r - 0.00295*cos(gm2 - em4) - 0.00263*cos(em + gs) + 0.00249*cos(gm3 - em2);
    r = r - 0.00221*cos(gm + em2 - gs) + 0.00185*cos(fm2 - em2) - 0.00161*cos(2 * (em - gs));
    r = r + 0.00147*cos(gm + fm2 - em2) - 0.00142*cos(em4) + 0.00139*cos(gm2 - em2 + gs);
    
    rmoon = 6378.14 * (r - 0.00118 * cos(gm - em4 + gs) - 0.00116 * cos(gm2 + em2) ... 
            - 0.0011 * cos(gm2 - gs));

    % Geocentric Equatorial Right Ascension and Declination
    a = sin(plon) * cos(eps_0) - tan(plat) * sin(eps_0);
    b = cos(plon);
       
    ran = atan2(a, b);
    dec = asin(sin(plat) * cos(eps_0) + cos(plat) * sin(eps_0) * sin(plon));

    % Moon's geocentric position vector (km)
    pos = zeros(3, 1);

    pos(1) = rmoon * cos(ran)*cos(dec);
    pos(2) = rmoon * sin(ran)*cos(dec);
    pos(3) = rmoon * sin(dec);
    
end