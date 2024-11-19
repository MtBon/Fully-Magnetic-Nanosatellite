function fa = fundamental_arguments_a(t)

    % Compute the Fundamental Luni-Solar and Planetary Arguments,
    % in radians, according to the IAU 2010 conventions. 
    %
    % Parameters 
    % ----------
    %   t: double
    %       Time expressed in TDB Julian centuries since J2000
    %
    % Returns 
    % -------
    %   fa: double(14, 1)
    %       Fundamental arguments vector, in radians, whose elements are respectively: 
    %        - Ma: Mean anomaly of the Moon 
    %        - Sa: Mean anomaly of the Sun
    %        - um: Mean longitude of the Moon minus mean longitude of the ascending node `F`
    %        - Ds: Mean elongation of the Moon from the Sun 
    %        - Om: Mean longitude of the Moon's ascending node
    %        - l_Me: Mercury's mean heliocentric longitude
    %        - l_Ve: Venus's mean heliocentric longitude
    %        - l_Ea: Earth's mean heliocentric longitude
    %        - l_Ma: Mars's mean heliocentric longitude
    %        - l_Ju: Jupiter's mean heliocentric longitude
    %        - l_Sa: Saturn's mean heliocentric longitude
    %        - l_Ur: Uranus's mean heliocentric longitude
    %        - l_Ne: Neptune's mean heliocentric longitude
    %        - pa: General precession in longitude
    %
    % References
    % ----------
    %   [1] Luzum, B. and Petit G. (2012), The IERS Conventions (2010).

    % Arcseconds in a full circle 
    ARCSECTURN = 1296000.0;

    % Arcseconds to radians conversion factor
    arcsec2rad = pi/648000;

    % Mean anomaly of the Moon
    Ma = 485868.249036 + t*(1717915923.2178 + t*(31.8792 ...
         + t*(0.051635 - 0.00024470*t)));
    Ma = mod(Ma, ARCSECTURN) * arcsec2rad;

    % Mean anomaly of the Sun
	Sa = 1287104.793048 + t*(129596581.0481 + t*(-0.5532 ...
         + t*(0.000136 - 0.00001149*t)));
    Sa = mod(Sa, ARCSECTURN) * arcsec2rad;

	% Mean argument of the latitude of the Moon. 
	um = 335779.526232 + t*(1739527262.8478 + t*(-12.7512 ...
        + t*(-0.001037 + 0.00000417*t)));
    um = mod(um, ARCSECTURN) * arcsec2rad;

	% Mean longitude of the ascending node of the Moon. 
	Om = 450160.398036 + t*(-6962890.5431 + t*(7.4722 ...
        + t*(0.007702 - 0.00005939*t)));
    Om = mod(Om, ARCSECTURN) * arcsec2rad;

    % Mean elongation of the Moon from the Sun. 
    Ds = 1072260.703692 + t*(1602961601.2090 + t*(-6.3706 ...
        + t*(0.006593 - 0.00003169*t)));
    Ds = mod(Ds, ARCSECTURN) * arcsec2rad;
    
    % Two pi
    pi2 = 2*pi;
    
    % General accumulated precession in longitude
    pa = mod(t*(0.024381750 + 0.00000538691*t), pi2);

    % Mean heliocentric longitude of Mercury 
    l_Me = mod(4.402608842 + 2608.7903141574*t, pi2);

    % Mean heliocentric longitude of Venus
    l_Ve = mod(3.176146697 + 1021.3285546211*t, pi2);
    
    % Mean heliocentric longitude of Earth
    l_Ea = mod(1.753470314 + 628.3075849991*t, pi2);

    % Mean heliocentric longitude of Mars
    l_Ma = mod(6.203480913 + 334.0612426700*t, pi2);

    % Mean heliocentric longitude of Jupiter
    l_Ju = mod(0.599546497 + 52.9690962641*t, pi2);

    % Mean heliocentric longitude of Saturn
    l_Sa = mod(0.874016757 + 21.3299104960*t, pi2);
    
    % Mean heliocentric longitude of Uranus
    l_Ur = mod(5.481293872 + 7.4781598567*t, pi2);
     
    % Mean heliocentric longitude of Neptune
    l_Ne = mod(5.311886287 + 3.8133035638*t, pi2);

    % Fundamental arguments vector: 
    fa = [Ma, Sa, um, Ds, Om, l_Me, l_Ve, l_Ea, l_Ma, ...
          l_Ju, l_Sa, l_Ur, l_Ne, pa];

end