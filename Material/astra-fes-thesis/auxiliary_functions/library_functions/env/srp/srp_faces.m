function [a_srp, t_srp, cth] = srp_faces(r_sun, r_sat, A, nu, N_faces, R_faces, ...
                                    S_faces, NS, c_spe, c_dif, m)

    % Compute the Solar Radiation Pressure (SRP) acceleration and Torque 
    % by modelling the spacecraft as a collection of plates.
    % 
    % Parameters 
    % ----------
    %   r_sun: double (3, 1)
    %       Sun position vector with respect to the central body [m]
    %   r_sat: double (3, 1)
    %       Spacecraft position vector with respect to the central body [m]
    %   A: double (3, 3)
    %       Direction Cosine Matrix (DCM) from the inertial to the
    %       spacecraft body frame. 
    %   nu: double
    %       Shadow function.
    %   N_faces: double (3, NS)
    %       Faces normal directions in the body frame.
    %   R_faces: double (3, NS)
    %       Faces center of pressure wrt the spacecraft's center of mass in
    %       the body frame in [m].
    %   S_faces: double (1, NS)
    %       Area of each face in [m²].
    %   NS: integer
    %       Total number of faces.
    %   c_spe: double 
    %       Specular reflection coefficient.
    %   c_dif: double 
    %       Diffusive reflection coefficient. 
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
    
    a_srp = zeros(3, 1); 
    t_srp = zeros(3, 1); 

    % Ensure all vectors are column vectors in [m]
    r_sun = reshape(r_sun, 3, 1);
    r_sat = reshape(r_sat, 3, 1);

    % Compute Spacecraft to Sun Vector in the Body Frame in [m]
    r_sat2sun = A*(r_sun - r_sat); 

    % Compute Faces to Sun vectors in the Body Frame in [m]
    R_f2s = r_sat2sun*ones(1, NS) - R_faces;

    % Normalise the vectors
    U_f2s = R_f2s; 
    for i = 1:NS 
        U_f2s(:, i) = U_f2s(:, i)/norm(U_f2s(:, i));
    end

    % Compute cosine of the incident angle! 
    cth = dot(U_f2s, N_faces);

    % If the cosine is negative, the face is not facing the Sun, 
    % thus it should not contribute to the acceleration
    cth = 0.5*(abs(cth) + cth);
    
    if nu == 0
        return; 
    end

    % Compute the SRP force on each face
    Fi_srp = zeros(3, NS);

    for i = 1:NS
        F_s = (1 - c_spe)*U_f2s(:, i); % Incident force
        F_n = 2*(c_dif/3 + c_spe*cth(i))*N_faces(:, i);  % Normal force

        % Total contribution
        Fi_srp(:, i) = S_faces(i)*cth(i)*(F_s + F_n);
    end

    % Compute solar pressure at given distance [Pa]
    d_sun = norm(r_sat2sun); % in [m]
    p = solar_pressure(1e-3*d_sun);
    
    % SRP Force [N] and Torque [Nm] on each face in the Body Frame 
    Fi_srp = -nu*p*Fi_srp;
    Ti_srp = cross(R_faces, Fi_srp);

    % Total Force [N] and Torque [Nm] in the Body Frame
    F_srp = sum(Fi_srp, 2);
    t_srp = sum(Ti_srp, 2);

    % SRP acceleration in the inertial frame [m/s²]
    a_srp = A'*F_srp/m;

end