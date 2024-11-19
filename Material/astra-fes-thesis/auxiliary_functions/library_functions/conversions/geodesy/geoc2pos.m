function pos = geoc2pos(r, lon, lat)
    
    % Transform geocentric coordinates in a cartesian position vector 
    % in the body-fixed frame. 
    %
    % Parameters
    % ----------
    %   r: double 
    %       Geocentric radial distance [km]
    %   lon: double 
    %       Geocentric longitude [rad]
    %   lat: double 
    %       Geocentric latitude [rad]
    % 
    % Returns  
    % -------
    %   pos: double (3, 1)
    %       Cartesian geocentric position vector. 

    slon = sin(lon); clon = cos(lon); 
    slat = sin(lat); clat = cos(lat); 

    pos = [r*clat*clon; r*clat*slon; r*slat];

end