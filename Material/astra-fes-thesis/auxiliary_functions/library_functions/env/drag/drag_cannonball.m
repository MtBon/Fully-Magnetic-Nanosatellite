function a_drag = drag_cannonball(pos, vel, CD, BC, rho)

    % Compute the Drag acceleration using the Cannonball approximation 
    % in [m/s²]. This function is valid only for the Earth.
    %
    % Parameters
    % ----------
    %   pos: double (3, 1)
    %       Spacecraft position vector in [m].
    %   vel: double (3, 1)
    %       Spacecraft velocity vector in [m/s]
    %   CD: double
    %       Drag Coefficient 
    %   BC: double 
    %       Ballistic Coefficient [kg/m²]
    %   rho: double 
    %       Air density [kg/m³]
    %
    % Returns
    % -------
    %   a_drag: double (3, 1)
    %       Drag acceleration in [m/s²]
    %
    % References
    % ----------
    % [1] H. D. Curtis, Orbital Mechanics for Engineering Students, IV Edition.
   

    % Ensure all vectors are column vectors and in [m]
    pos = reshape(pos, 3, 1);
    vel = reshape(vel, 3, 1);

    % Earth's angular velocity [rad/s]
    we = 7.292115486e-5;

    % Relative velocity wrt the atmosphere
    vrel = vel - cross3([0; 0; we], pos);
    
    % Drag acceleration [m/s]
    a_drag = - 0.5*rho*CD/BC*norm(vrel)*vrel;

end