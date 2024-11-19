function nu = shadow_function(r_sat, r_sun, Rs, Rp)

    % Compute the illumination conditions of a spacecraft given its
    % position vector, the sun position and radius and the radius of the
    % occulting body. This function can only be used when the occulting
    % body is smaller than the sun. When this assumption is not verified,
    % some issues could arise because umbra zones between the sun and the
    % body are not checked. The occulting body is assumed to be spherical.
    %
    % Parameters
    % ----------
    %  r_sat: double(3, 1)
    %       Position vector of the spacecraft wrt the occulting body.
    %  r_sun: double(3, 1)
    %       Position vector of the Sun wrt the occulting body.
    %  Rs: double 
    %       Radius of the Sun 
    %  Rp: double 
    %       Radius of the occulting body 
    %
    % Returns
    % -------
    %  nu: integer 
    %       Dimensionless shadow function coefficient, between 0 and 1. 
    %       Which expresses the luminosity ratio during eclipses.    
    %
    % Notes 
    % -----
    %   All the inputs must be expressed in the same unit of measure.
    %
    % References
    % ----------
    % [1] O. Montenbruck, E. Gill, Satellite Orbits Models, Methods and
    %     Applications, 2005, Springer

    % Ensure both vectors are column vectors 
    r_sat = reshape(r_sat, 3, 1); 
    r_sun = reshape(r_sun, 3, 1); 
    
    % Magnitude of the spacecraft and sun vectors 
    n_sat = sqrt(r_sat(1)^2 + r_sat(2)^2 + r_sat(3)^2); 
    n_sun = sqrt(r_sun(1)^2 + r_sun(2)^2 + r_sun(3)^2);

    % Retrieve unit vectors 
    u_sat = r_sat/n_sat;
    u_sun = r_sun/n_sun;

    % Cross product w = (u_sat, u_sun): 
    w1 = u_sat(2)*u_sun(3) - u_sat(3)*u_sun(2);
    w2 = u_sat(3)*u_sun(1) - u_sat(1)*u_sun(3);
    w3 = u_sat(1)*u_sun(2) - u_sat(2)*u_sun(1);

    % Retrieve the sin and cos of the angle between u_sat and u_sun
    sk = sqrt(w1*w1 + w2*w2 + w3*w3);
    ck = u_sat(1)*u_sun(1) + u_sat(2)*u_sun(2) + u_sat(3)*u_sun(3);
    
    % Find vertical and horizontal projections of the spacecraft
    % position onto the shadow axis 
    sat_v = n_sat*sk; 
    sat_h = n_sat*abs(ck); % ck could be negative !
    
    % Compute the penumbra angle
    sfp = (Rs + Rp) / n_sun;

    fp = asin(sfp);
    x = Rp/sfp;

    if ck > 0 
        % The spacecraft is between the sun and the planet 
        % Perform a check for a possible penumbra zone! 
        
        % Penumbra vertical size at spacecraft location
        pen_v = tan(fp)*(x - sat_h);

        % Horizontal point at which penumbra starts 
        pen_x = Rp*sfp;

        if sat_v > pen_v || sat_h > pen_x 
            nu = 1; % No eclipse! 
        else 
            nu = compute_eclipse(r_sat, n_sat, u_sat, r_sun, Rs, Rp);
        end 
        
        return; 
    end

    % The spacecraft is beyond the planet

    % Penumbra vertical size at spacecraft location 
    pen_v = tan(fp)*(x + sat_h);

    if sat_v > pen_v 
        nu = 1; % No eclipse!
        return;
    end
    
    % At this point the satellite is either in penumbra or umbra

    % Compute the umbra angle 
    sfu = (Rs - Rp) / n_sun; 

    % Check for the penumbra portion at the first mid-half 
    if sat_h < Rp*sfu 
        nu = compute_eclipse(r_sat, n_sat, u_sat, r_sun, Rs, Rp);
        return;
    end

    % Retrieve the umbra angle
    fu = asin(sfu);
    y = Rp/sfu;

    % Compute umbra vertical size at spacecraft location
    umb_v = tan(fu)*abs(sat_h - y);

    if sat_v < umb_v && sat_h < y 
        nu = 0; % Total eclipse
    else
        % Penumbra or Annular eclipse
        nu = compute_eclipse(r_sat, n_sat, u_sat, r_sun, Rs, Rp);
    end

end

function nu = compute_eclipse(r_sat, n_sat, u_sat, r_sun, Rs, Rp)

    % Compute the luminosity ratio between eclipses 
    % See Montenbruck for references
    
    % Compute relative sun-sc position vector 
    r_rel = r_sun - r_sat; 
    n_rel = sqrt(r_rel(1)^2 + r_rel(2)^2 + r_rel(3)^2);
    
    u_rel = r_rel/n_rel; 
    
    % Compute apparent sun and planet radiuses as seen by the spacecraft
    a = asin(Rs/n_rel); 
    b = asin(Rp/n_sat);

    % Compute the angular separation between the planet and the sun 
    c = acos(dot(-u_sat, u_rel));
    
    if abs(a-b) < c && c < a + b
        % Partial eclipse
        a2 = a^2; 
        b2 = b^2; 
        
        x = (c^2 + a2 - b2)/(2*c);
        y = sqrt(a2 - x^2); 
        A = a2*acos(x/a) + b2*acos((c-x)/b) - c*y; 
        
        nu = 1 - A/(pi*a2);
        
    else 
        % Annular eclipse
        nu = 1 - (1-cos(b))/(1-cos(a));
        
    end
end