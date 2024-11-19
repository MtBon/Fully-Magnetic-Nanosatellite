function [a_drag, t_drag] = drag_faces(pos, vel, w, A, rho, N_faces, R_faces, ... 
                                S_faces, NS, CD, m)

    % Compute the drag acceleration and Torque by modelling the spacecraft 
    % as a collection of plates. It is valid only for the Earth.
    % 
    % Parameters 
    % ----------
    %   pos: double (3, 1)
    %       Spacecraft position vector with respect to the central body [m]
    %   vel: double (3, 1)
    %       Spacecraft velocity vector with respect to the central body [m/s]
    %   w: double (3, 1)
    %       Spacecraft angular velocity vector in the Body Frame [rad/s]
    %   A: double (3, 3)
    %       Direction Cosine Matrix (DCM) from the inertial to the
    %       spacecraft body frame. 
    %   rho: double
    %       Air density [kg/m³]
    %   N_faces: double (3, NS)
    %       Faces normal directions in the body frame.
    %   R_faces: double (3, NS)
    %       Faces center of pressure wrt the spacecraft's center of mass in
    %       the body frame in [m].
    %   S_faces: double (1, NS)
    %       Area of each face in [m²].
    %   NS: integer
    %       Total number of faces.
    %   CD: double 
    %       Drag coefficient.
    %   m: double 
    %       Total spacecraft mass [kg].
    %
    % Returns 
    % -------
    %   a_srp: double (3, 1)
    %       Total SRP acceleration in the inertial frame in [m/s²].
    %   t_srp: double (3, 1)
    %       Total SRP torque in the body frame in [Nm].
    %   cth: double (1, NS)
    %       Cosines of the incident angles between each face normal and 
    %       the Sun's direction. Set to 0 if in shadow. 
    %
    % References
    % ----------
    %   [1] D. Vallado, Fundamentals of Astrodynamics, IV Edition, p. 582

    % Ensure all vectors are column vectors
    pos = reshape(pos, 3, 1);
    vel = reshape(vel, 3, 1);
    w   = reshape(w, 3, 1);

    % Earth's angular velocity [rad/s]
    we = 7.292115486e-5;

    % SC Relative velocity wrt the atmosphere in the Body Frame [m/s]
    vrel = A*(vel - cross([0; 0; we], pos));

    % Relative faces velocites wrt the atmosphere in the Body Frame [m/s]
    v_faces = vrel*ones(1, NS) + cross(w*ones(1, NS), R_faces);

    % Normalise each incident direction
    u_faces = v_faces; 
    for i = 1:NS 
        u_faces(:, i) = u_faces(:, i)/norm(u_faces(:, i));
    end

    % Find which faces are exposed to the air flow
    % i.e., the cosine of the incident angle > 0 
    cth = dot(u_faces, N_faces);
    cth = 0.5*(abs(cth) + cth);

    % Compute the drag force acting on each face in the Body Frame
    Fi_drag = zeros(3, NS); 
    for i = 1:NS 
        % cth accounts for the face area projection along the relative velocity direction.
        Fi_drag(:, i) = S_faces(i)*cth(i)*norm(v_faces(:, i))*v_faces(:, i);
    end

    % Drag Force and Torque on each face in the Body Frame
    Fi_drag = -0.5*rho*CD*Fi_drag;
    Ti_drag = cross(R_faces, Fi_drag); 

    % Total Force [N] and Torque [Nm] in the Body Frame
    F_drag = sum(Fi_drag, 2);
    t_drag = sum(Ti_drag, 2);

    % Drag acceleration in the inertial frame [m/s²]
    a_drag = A'*F_drag/m;

end