function dat = get_leapseconds(utc, leap)

    % For a given UTC date, extract the number of leapseconds
    %
    % Parameters
    % ----------
    %   utc: number
    %       UTC Julian Days since J2000.0
    %   leap: number(M, 2)
    %       Leap seconds look-up table, computed from the function 
    %       load_leapseconds. 
    %
    % Returns 
    % -------
    %   dat: number
    %       TAI - UTC offset, in seconds.
    
    %#codegen
    
    dat = leap(find(leap(:, 1) <= utc, 1, 'last'), 2);
 
end


