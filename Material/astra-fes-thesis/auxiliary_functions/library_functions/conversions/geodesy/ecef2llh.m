function [lat, lon, h] = ecef2llh(pos)

    % Compute geodetic latitude, longitude and altitude from the reference
    % WGS-84 ellipsoid. 
    %
    % Parameters
    % ----------
    %   pos: double (3, 1)
    %       Position vector in the ECEF frame in [m]. 
    % 
    % Returns 
    % -------
    %   lat: double 
    %       Geodetic latitude [rad]. 
    %   lon: double 
    %       Geocentric longitude [rad]. 
    %   h: double 
    %       Altitude over the reference ellipsoid [m]. 

    % Reference WGS-84 ellipsoid parameters

    % Equatorial radius [m];
    R = 6378137; 

    % Ellipsoid Flattening 
    f = 1 / 298.257;
    
    [lat, lon, h] = pos2geod(pos, R, f);

end