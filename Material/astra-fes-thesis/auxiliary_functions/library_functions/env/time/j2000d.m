function d = j2000d(jd)

    % Convert a Julian Date (in days) to days since J2000.0
    % 
    % Parameters
    % ----------
    %   jd: double
    %       Julian Date [days]
    %  
    % Returns 
    % -------
    %   d:  double
    %       Julian Days since the reference epoch J2000.0
    
    d = jd - 2451545.0;

end