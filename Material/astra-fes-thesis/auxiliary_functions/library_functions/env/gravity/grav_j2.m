function a_grav = grav_j2(r_sat, Rb, GM, j2)

    % Compute the gravitational J2 perturbation at a given position. 
    %
    % Parameters
    % ----------
    %   r_sat: double (3, 1)
    %       Satellite position vector 
    %   Rb: double 
    %       Central body radius
    %   GM: double 
    %       Planetary gravitational constant.
    %   j2: J2 coefficient
    %
    % Returns
    % -------
    %   a_grav: double (3, 1)
    %       J2 gravitational perturbation. 
    %
    % Notes
    % -----
    %   The input units of length must be coherent (i.e., if r_sat is
    %   expressed in km, then Rb and GM will aswell). The output
    %   acceleration is expressed in the same length unit of the inputs.

    r2 = r_sat(1)^2 + r_sat(2)^2 + r_sat(3)^2;
    r = sqrt(r2);
    r3 = r*r2;
    r5 = r3*r2;

    a = 1.5*j2*GM*Rb^2/r5;
    b = 5*(r_sat(3)/r)^2;
    
    a_grav = [a*(b - 1)*r_sat(1); 
              a*(b - 1)*r_sat(2);
              a*(b - 3)*r_sat(3)];
    
end
