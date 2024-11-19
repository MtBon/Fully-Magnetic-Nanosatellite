function mjd = jd2mjd(jd)

    % Convert Julian Days to Modified Julian Days
    %
    % Parameters
    % ----------
    %   jd: double 
    %       Julian Date [days]
    %   
    % Returns 
    % -------
    %   mjd: double
    %       Modified Julian Date [days]

    mjd = jd - 2400000.5;

end