   
function [g, h] = get_magnetic_coeffs(time, mag_data)

    % Compute IGRF13's magnetic field coefficients at epoch 
    % 
    % Parameters
    % ----------
    %   time: double
    %       Terrestrial Time (TT) seconds since J2000.0. 
    %   mag_data: double (maxdeg, 6)
    %       Matrix storing the IGRF model orders and magnetic coefficients.
    %       Computed with teh function `load_mag_model`
    %
    % Returns 
    % -------
    %   g: double (maxdeg, maxdeg+1)
    %       Gaussian normalised cosine coefficients at epoch.
    %   h: double (maxdeg, maxdeg+1)
    %       Gaussian normalised sine coefficients at epoch. 
    
    
    % It doesn't really matter which timescale we are using here, since the
    % differences between them are in the order of seconds, whereas the
    % IGRF coefficients have yearly variations. 
    
    % Reference epoch deve essere 01/01/2020 00:00:00 se si usa IGRF20!
    % Compute days and years since 2020 at midnight
    days = time/86400 - 7304.5;
    t_yrs = days/365.25;
    
    % Maximum IGRF13 expansion degree
    maxdeg = 13;

    g = zeros(maxdeg, maxdeg + 1); 
    h = zeros(maxdeg, maxdeg + 1);
    
    % Compute Gaussian normalised magnetic coefficients at current epoch!
    for n = 1:maxdeg 
        for m = 0:n
            idx = n*(n+1)/2 + m;
        
            g(n, m+1) = mag_data(idx, 3) + mag_data(idx, 4)*t_yrs; 
            h(n, m+1) = mag_data(idx, 5) + mag_data(idx, 6)*t_yrs; 
        end
    end
    
end