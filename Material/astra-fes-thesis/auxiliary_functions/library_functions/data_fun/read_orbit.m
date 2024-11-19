function [fid, oev] = read_orbit(filename)

    % Read Classical Orbital Elements from a datafile named `filename`. 
    % 
    % Inputs 
    % ------
    %   filename: string
    %       Name of the data file containing the orbital elements. It is
    %       assumed that the angular measurements are expressed in degrees.
    %  
    % Outputs
    % -------
    %   fid: integer
    %       File ID
    %   oev: double(6, 1)
    %       Classical Orbital Elements (COE) vector. All the angular
    %       elements are expressed in radians. The semi-major axis is
    %       reported in the same units of the data-file. The COE order 
    %       is: [sma, ecc, inc, ran, aop, tan]
    
    % Degrees to Radians converter
    dtr = pi / 180;
    
    % Open data file
    fid = fopen(filename, 'r');
    
    % Check for opening error
    while fid == -1
        clc
        
        fprintf('\n\nError: cannot find the orbital file!\n\n');
        [filename, path] = uigetfile('*.in', 'please select an input data file');
        fid = fopen(fullfile(path, filename), 'r');

    end
    
    % read 23 lines of data file
    for i = 1:1:23
        
        cline = fgetl(fid);
        
        switch i
            case 3 % Semi-major axis
                oev(1) = str2double(cline);
                
            case 7 % Eccentricity
                oev(2) = str2double(cline);
                
            case 11 % Inclination
                oev(3) = dtr * str2double(cline);
                
            case 15 % Argument of Perigee 
                oev(5) = dtr * str2double(cline);
                
            case 19 % Right-Ascension of the Ascending Node
                oev(4) = dtr * str2double(cline);
                
            case 23 % True Anomaly
                oev(6) = dtr * str2double(cline);
                
        end
    end
    
    status = fclose(fid); %#ok<NASGU>

end
