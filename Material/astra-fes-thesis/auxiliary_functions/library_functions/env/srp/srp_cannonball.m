function a_srp = srp_cannonball(r_sat, r_sun, Cr, BC, nu)

    % Compute the Solar Radiation Pressure acceleration using 
    % the Cannonball approximation in [m/s²]. It is valid for any 
    % central body.
    %
    % Parameters
    % ----------
    %   r_sat: double (3, 1)
    %       Spacecraft position vector with respect to the central body [m].
    %   r_sun: double (3, 1)
    %       Sun position vector with respect to the central body [m]
    %   Cr: double
    %       Reflectivity coefficient.
    %   BC: double 
    %       Ballistic Coefficient [kg/m²]
    %   nu: double 
    %       Shadow function.
    %
    % Returns
    % -------
    %   a_srp: double (3, 1)
    %       SRP acceleration in [m/s²]
    %
    % References
    % ----------
    % [1] H. D. Curtis, Orbital Mechanics for Engineering Students, IV Edition.
   
    
    r_sat = reshape(r_sat, 3, 1); 
    r_sun = reshape(r_sun, 3, 1); 
    
    % Compute Sun to Satellite relative vector
    dr = r_sat - r_sun;
    d_sun = norm(dr);
    
    u = dr/d_sun;
    
    % Compute solar pressure at given distance [Pa]
    p = solar_pressure(1e-3*d_sun);
    
    % Compute Solar Radiation Pressure acceleration 
    a_srp = nu*p*Cr/BC*u; % [m/s²] 
    
end