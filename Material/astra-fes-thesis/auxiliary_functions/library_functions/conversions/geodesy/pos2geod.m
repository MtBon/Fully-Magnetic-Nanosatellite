function [lat, lon, h] = pos2geod(pos, R, f)

    % Transform a cartesian body-fixed position vector into longitude,
    % geodetic latitude and altitude over the reference ellipsoid.
    %
    % Parameters
    % ----------
    %   pos: double (3, 1)
    %       Position vector in the body-fixed frame [km].
    %   R: double 
    %       Reference body radius [km].
    %   f: double 
    %       Body flattening. 
    %   
    % Returns 
    % -------
    %   lat: double 
    %       Geodetic latitude [rad]. 
    %   lon: double 
    %       Geocentric longitude [rad]. 
    %   h: double 
    %       Altitude over the reference ellipsoid [km]. 
    %
    % References 
    % ----------
    %   [1] D. A. Vallado, Fundamentals of Astrodynamics and Applications

    % Iterative algorithm settings
    tol = 1e-11;
    maxiter = 5;
    
    x = pos(1); 
    y = pos(2); 
    z = pos(3);

    sz = sign(z); 
    
    % Distance from the polar axis
    rho = x^2 + y^2; 

    % Square of the eccentricity
    e2 = f*(2-f); 

    zn = z; 
    cn = 1.0; 
    sphi2 = 0.0;

    % Iteratively solve for the geodetic latitude. Note the equations are
    % efficiently written to avoid any trigonometric function. 
    
    iter = 0; 
    err = tol + 1; 
    while err > tol && iter < maxiter 
        c = cn; 
        z2 = zn^2; 
        sphi2 = z2/(rho + z2);
        
        cn = R*e2*sqrt(sphi2/(1-e2*sphi2));
        zn = z + cn*sz; 

        err = abs(cn - c);

        iter = iter + 1; 
    end

    lon = atan2(y, x); 
    lat = atan(zn/sqrt(rho));
    h = sqrt(rho + zn^2) - R/sqrt(1-e2*sphi2);

end