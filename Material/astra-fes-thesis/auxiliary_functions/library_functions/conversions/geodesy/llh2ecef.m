function pos = llh2ecef(lat, lon, h)

    % Transform geodetic latitude, longitude and altitude from the WGS-84
    % ellipsoid into a position vector in the ECEF frame in [km].
    %
    % Parameters 
    % ----------
    %   lat: double 
    %       Geodetic latitude [rad]. 
    %   lon: double 
    %       Geodetic longitude [rad].
    %   h: double
    %       Altitude from the reference WGS-84 ellipsoid [m]. 
    %   
    % Returns
    % -------
    %   pos: double (3, 1)
    %       Position vector in the ECEF frame in [m].

    % Equatorial radius [m];
    R = 6378137; 

    % Ellipsoid Flattening 
    f = 1 / 298.257;

    pos = geod2pos(lat, lon, h, R, f);

end

