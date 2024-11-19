function jd = mjd2jd(mjd)

    % Convert Modified Julian Days to Julian Days
    %
    % Parameters
    % ----------
    %   mjd: double 
    %       Modified Julian Date [days]
    %   
    % Returns 
    % -------
    %   jd: double
    %       Julian Date [days]

    jd = mjd + 2400000.5;

end