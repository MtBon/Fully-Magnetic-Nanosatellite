function tle = read_tle_lines(filename_tle)

    % Read a TLE file and return its lines
    %
    % Parameters
    % ----------
    %   filename_tle: string
    %       TLE filename
    %
    % Returns
    % -------
    %   tle: struct
    %       TLE structure containing the three TLE lines.
    
    tleID = fopen(filename_tle, 'r');
    tle.line0 = string(fgetl(tleID));
    tle.line1 = string(fgetl(tleID));
    tle.line2 = string(fgetl(tleID));

end
