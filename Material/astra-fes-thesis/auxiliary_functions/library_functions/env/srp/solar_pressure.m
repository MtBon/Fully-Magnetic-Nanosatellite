function p = solar_pressure(d_sun)
    
    % Compute the solar radiation pressure in W/m² at the 
    % given distance from the sun. 
    %
    % Parameter
    % ---------
    %   d_sun: double 
    %       Distance from the sun in [km].
    %
    % Returns
    % -------
    %   p: double 
    %       Solar radiation pressure in [Pa]
    %
    % References
    % ----------
    %   [1] H. D. Curtis, Orbital Mechanics for Engineering Students, IV Edition.
    
    % Photosphere radius [km] 
    R0 = 696000; 
    
    % Radiated power intensity at the sun's surface [W/m²]
    S0 = 63.15e6; 
    
    % Speed of light [m/s]
    c = 299792458; 
    
    % Compute the energy flux [W/m²]: 
    S = S0 * (R0/d_sun)^2;
    
    % Compute the solar radiation pressure: 
    p = S/c;

end