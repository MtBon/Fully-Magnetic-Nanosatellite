function LEAP = load_leapseconds(filename)
    
    % Extract the LEAP seconds array from the NAIF Leapseconds kernel. 
    %
    % Parameters
    % ----------
    %   filename: str 
    %       Name of the leapsecond kernel file 
    %   
    % Returns 
    % -------
    %   LEAP: double (n, 2)
    %       Leap seconds array. The first column is the JD days since J2000 at
    %       which the leapseconds (in the second column) have been updated.
    
    % Pre-allocated a large enough array
    LEAP = zeros(100, 2);

    % Julian days at 01-01-2000 12:00:00
    DJ2000 = 2451545;

    % Regex expression to extract leap seconds rows
    re = "(?<dat>[0-9]{2}),\s+@(?<date>[0-9]{4}-[A-Z]{3}-[0-9])";

    fileID = fopen(filename, 'r');

    nleap = 0;
    line = fgetl(fileID);
    while ischar(line)

        % Extract only the leap seconds data
        idx = regexp(line, re, 'ONCE');
        if ~isempty(idx)
            % Match regex expression! 
            tmp = regexp(line, re, 'match');
            [data, ~] = regexp(tmp, ',\s+@', 'split');

            nleap = nleap + 1; 

            % Save the leapsecond date as Julian Days since J2000.0
            dtime = datetime(data{1}{2}, 'InputFormat', 'yyyy-MMM-dd'); 
            LEAP(nleap, 1) = juliandate(dtime) - DJ2000; 

            % Leap second value
            LEAP(nleap, 2) = str2double(data{1}{1});
        end
        
        % Proceed to the next line
        line = fgetl(fileID);
    end
    
    % Clear all the preassigned empty rows of LEAP
    LEAP(nleap+1:end, :) = [];
    fclose(fileID);
    
end



