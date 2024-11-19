function [sim_start, sim_duration] = read_time_settings(filename_time)

    % Retrieve the simulation time parameters
    %
    % Parameters
    % ----------
    %   filename_time: string
    %       Name of the file with the simulation time settings.
    %
    % Returns 
    % -------
    %   sim_start: double
    %       Julian Date at simulation start time 
    %   sim_duration: double
    %       Simulation duration, in seconds
    
    fid = fopen(filename_time, 'r');
    
    % Check for file open error
    while fid == -1
        clc
        
        fprintf('\n\nError: cannot find the parameter file!\n\n');
        [filename, path] = uigetfile('*.in', 'please select an input data file');
        fid = fopen(fullfile(path,filename), 'r');
        
    end

    % Read 29 lines of data file
    for i = 1:1:19
        cline = fgetl(fid);
        switch i
            case 7
                month = str2double(cline);
            case 9
                day = str2double(cline);
            case 11
                year = str2double(cline);
            case 15
                hour = str2double(cline);
            case 17
                min = str2double(cline);
            case 19
                sec = str2double(cline);
            case 2
                ndays = str2double(cline);
        end
    end
    
    % Retrieve starting Julian Date
    sim_start = juliandate(year, month, day, hour, min, sec);

    % Compute simulation duration
    sim_duration = 86400*ndays;

end