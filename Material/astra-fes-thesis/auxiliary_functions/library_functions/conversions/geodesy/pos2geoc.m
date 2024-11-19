function [lat, lon, r] = pos2geoc(pos)

    % Convert a cartesian position vector in the fixed-body frame into geocentric
    % latitude, longitude and radius.
    %
    % Parameters
    % ----------
    %   pos: double (3, 1)
    %       Geocentric position vector [km]
    %  
    % Returns 
    % -------
    %   lat: double 
    %       Geocentric latitude [rad].
    %   lon: double 
    %       Geocentric longitude [rad].
    %   r: double 
    %       Geocentric radial distance [km]
    

    r = sqrt(pos(1)^2 + pos(2)^2 + pos(3)^2);

    lon = atan2(pos(2), pos(1)); 
    lat = asin(pos(3)/r); 

end