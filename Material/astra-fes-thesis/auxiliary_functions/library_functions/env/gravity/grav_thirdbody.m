function acc = grav_thirdbody(GM, r_sc, r_body)

    % Compute the third body gravitational perturbation in an 
    % inertial reference frame. This function exploits a different formulation 
    % to avoid numerical problems when the spacecraft is very far from the 
    % given body. The result is unchanged.
    % 
    % Parameters
    % ----------
    %   GM: double 
    %        Third body gravitational constant.
    %   r_sc: double (3, 1)
    %        Spacecraft inertial position vector.
    %   r_body: double (3, 1)
    %        Third body inertial position vector.
    %
    % Notes
    % -----
    %   r_sc and r_body must be expressed in the same inertial reference frame. 
    %
    % References
    % ----------
    %   [1] H. D. Curtis, Orbital Mechanics for Engineering Students, IV Edition.

    q = dot(r_sc, 2*r_body - r_sc)/dot(r_body, r_body);
    f = q*(q^2 - 3*q + 3)/(1 + (1 - q)^1.5);

    dr = r_sc - r_body; 

    acc = GM*(f*r_body - r_sc)/norm(dr)^3;

end