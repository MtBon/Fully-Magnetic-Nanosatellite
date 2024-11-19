function coe = car2coe(stv, GM)

    % Convert a cartesian state (position, velocity) into Classical 
    % Orbital Elements (sma, ecc, inc, ran, aop, tan). 
    % 
    % Parameters
    % ----------
    %
    % Returns 
    % -------
    %
    %

    %#codegen

    tol = 1e-10; 
    pos = stv(1:3, :); 
    vel = stv(4:6, :); 

    % Compute angular momentum
    h = cross(pos, vel);
    hn = norm(h); 

    % Compute line of nodes
    N = cross([0; 0; 1], h);

    % Comput eccentricity vector  
    ecc_vec = eccentricity_vector(stv, GM);

    % Compute eccentricity 
    ecc = norm(ecc_vec); 

    % Compute semi-major axis 
    p = hn.^2/GM;
    sma = p./(1-ecc.^2); 

    % Compute inclination 
    inc = acos(h(3, :)./hn);

    % Check which orbits are circular and/or equatorial
    iscircular = ecc < tol; 
    isequatorial = inc < tol; 

    aop = zeros(size(sma));
    ran = zeros(size(sma));
    theta = zeros(size(sma)); 

    % Handle equatorial non-circular orbits
    i = isequatorial & ~iscircular;
    if any(i)

        % In this case we compute the true longitude of periapsis, 
        % which combines ran and aop to remove the ambiguity
        aop(i) = atan2(ecc_vec(2, i), ecc_vec(1, i)); 
        
        % These values are both not normalised by ecc_vec and pos norms
        sth = dot(h(:, i), cross(ecc_vec(:, i), pos(:, i)))./hn(i);
        cth = dot(pos(:, i), ecc_vec(:, i));

        theta(i) = atan2(sth, cth);

    end

    % Handle circular inclined orbits
    i(:) = ~isequatorial & iscircular;
    if any(i)

        ran(i) = atan2(N(2, i), N(1, i)); 

        % For a circular orbit, we use in-place of the 
        % true anomaly, the argument of latitude, i.e., ω + θ
        
        % These are both not normalised by pos and N, so there is no issue
        su = dot(pos(:, i), cross(h(:, i), N(:, i)))./hn(i);
        cu = dot(pos(:, i), N(:, i)); 

        theta(i) = atan2(su, cu);

    end

    % Handle circular equatorial orbits 
    i(:) = iscircular & isequatorial;
    if any(i)
        % For these orbits, we define theta = true longitude, 
        % i.e., the angle measure eastward from the x-axis to 
        % the position of the satellite
        theta(i) = atan2(pos(2, i), pos(1, i));
    end

    % Handle all remaining orbits
    i(:) = ~iscircular & ~isequatorial;
    if any(i)

        ran(i) = atan2(N(2, i), N(1, i)); 

        % Compute true anomaly
        r = norm(pos(:, i));
        theta(i) = acos(dot(ecc_vec(:, i)./ecc(i), pos(:, i)./r));
        
        % Compute radial velocity 
        vr = dot(pos, vel); 
        
        % Quadrant check! 
        idx = i & vr < 0; 
        theta(idx) = 2*pi - theta(idx); 

        % Argument of perigee is computed from argument of latitude 
        slat = dot(pos(:, i), cross(h(:, i), N(:, i)))./hn(i);
        clat = dot(pos(:, i), N(:, i)); 

        aop(i) = atan2(slat, clat) - theta(i); 

    end

    coe = [sma; ecc; inc; ...
           mod(ran,2*pi); ...
           mod(aop,2*pi); ...
           mod(theta,2*pi)];
end