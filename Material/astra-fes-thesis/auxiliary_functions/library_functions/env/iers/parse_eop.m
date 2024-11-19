function EOP = parse_eop(filename, jdstart, jdend)
    
    if jdstart > jdend
        error("Start Julian Date greater than ending date.");
    end    
    
    % Transform input julian days to Modified Julian Days
    mjd_start = jd2mjd(jdstart);
    mjd_end = jd2mjd(jdend);
    
    opts = detectImportOptions(filename);
    eop_tab = readtable(filename, opts); 

    ncols = size(eop_tab, 1); 
    
    % Days Threshold to extract data ndays before and after the actual 
    % simulation start and end dates. 
    dthr = 5;
    
    col_start = floor(mjd_start) - table2array(eop_tab(1, 1)) + 1;
    col_end = ceil(mjd_end) - table2array(eop_tab(1, 1)) + 1;
    
    if col_start < 0 || col_end > ncols
        warning('EOP data are not valid for this date interval')
    end
    
    idx_start = max(1, col_start - dthr);
    idx_stop = min(ncols, col_end + dthr);
    
    % Desired EOP Data columns to be extracted
    idx = idx_start:idx_stop;
    
    % Arcseconds to radians conversion factor
    a2r = arcsec2rad(1);
    
    % Extract and convert EOP data to [s] and [rad]
    EOP.dut1 = table2array(eop_tab(idx, 11));               % in [s]
    EOP.xp = a2r*table2array(eop_tab(idx, 6));              % in [rad]
    EOP.yp = a2r*table2array(eop_tab(idx, 7));              % in [rad]
    EOP.LOD = 1e-3*table2array(eop_tab(idx, 13));           % in [s]
    EOP.dx = 1e-3*a2r*table2array(eop_tab(idx, 20));        % in [rad]
    EOP.dy = 1e-3*a2r*table2array(eop_tab(idx, 21));        % in [rad]
    
    EOP.jd = j2000d(mjd2jd(table2array(eop_tab(idx, 1))));  % in J2000 [days]
    
end
