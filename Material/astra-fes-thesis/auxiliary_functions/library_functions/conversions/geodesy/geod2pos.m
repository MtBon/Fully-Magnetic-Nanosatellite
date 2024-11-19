function pos = geod2pos(lat, lon, h, R, f)

    % Transform geodetic longitude, latitude and altitude from the
    % reference ellipsoid into a cartesian body-fixed position vector. 
    %
    % Parameters
    % ----------
    %   lat: double
    %       Geodetic latitude [rad].
    %   lon: double 
    %       Geodetic longitude [rad]. 
    %   h: double
    %       Altitude from the reference ellipsoid [km]. 
    %   R: double
    %       Reference body radius [km]. 
    %   f: double
    %       Flattening parameter.
    %   
    % Returns
    % -------
    %   pos: double (3, 1)
    %       Cartesian body-fixed position vector. 
    %
    % References
    % ----------
    %   [1] D. A. Vallado, Fundamentals of Astrodynamics and Applications

    % Compute square of eccentricity from flattening
    e2 = (2-f)*f; 
    
    slat = sin(lat); clat = cos(lat);
    slon = sin(lon); clon = cos(lon); 

    d = R/sqrt(1-e2*slat^2);
    c = (d + h)*clat; 
    s = (1 - e2)*d; 

    pos = [c*clon; c*slon; (s+h)*slat];

end