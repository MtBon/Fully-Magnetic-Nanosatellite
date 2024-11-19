function B = magnetic_field_spherical(pos, lat, lon, g, h)

    % Compute the magnetic field vector in the inertial frame using the 
    % spherical harmonic expansion and IGRF13's coefficients. 
    %
    % Parameters
    % ----------
    %   pos: double (3, 1)
    %       Spacecraft inertial position wrt the Earth [m].
    %   lat: double 
    %       Spacecraft geocentric latitude [rad].
    %   lon: double 
    %       Spacecraft East longitude from Greenwitch [rad].
    %
    % Returns
    % -------
    %   B: double (3, 1)
    %       Magnetic field vector in the inertial frame [T]. 
    %   
    % References
    % ----------
    %   [1]: J. R. Wertz, Spacecraft Attitude Determination and Control


    % Reference Earth geomagnetic radius [m]
    a = 6371200;
    
    % Radial distance [m]
    r = norm(pos);

    % Compute declination and right ascension!
    dec = asin(pos(3)/r);
    ras = atan3(pos(2), pos(1));

    maxdeg = size(g, 1);

    % Pre-compute the sin e cos terms to avoid repeated computations 
    [cth, sth] = sincos(maxdeg, lon);

    % Find co-elevation from geocentric latitude 
    th = pi/2 - lat;

    % Handle poles singularities!
    if abs(th) < 1e-8 
        th = 1e-8;
    elseif abs(th - pi) < 1e-8
        th = pi - 1e-8;
    end

    % Pre-compute the legendre polynomials 
    [P, dP] = compute_legendre(maxdeg, th);

    % Magnetic field along the radial, tangential and normal directions
    Br = 0.0; Bt = 0.0; Bp = 0.0; 

    % Ratio between a and current radial distance
    r_ratio = (a/r);
    r_ratio_n = r_ratio^2;

    for n = 1:maxdeg

        % Single components contributions at n-th deg.
        dbdr = 0.0; dbdt = 0.0; dbdp = 0.0;

        % Current order index
        crt = n*(n+1)/2 + 1;

        for m = 0:n 
            % Current index for legendre vector!
            idx = crt + m;

            % Current deg, ord magnetic coefficients.
            gnm = g(n, m+1); 
            hnm = h(n, m+1); 

            Pnm = P(idx);
            dPnm = dP(idx);
            
            tmp = gnm*cth(m+1) + hnm*sth(m+1);
            
            dbdr = dbdr + tmp*Pnm; 
            dbdt = dbdt + tmp*dPnm;
            dbdp = dbdp + m*(-gnm*sth(m+1) + hnm*cth(m+1))*Pnm;

        end

        % Update (a/r) ratio
        r_ratio_n = r_ratio_n*r_ratio; 
        
        % Finalise all the orders contributions of the n-th degree.
        Br = Br + (n+1)*r_ratio_n*dbdr;
        Bt = Bt + r_ratio_n*dbdt; 
        Bp = Bp + r_ratio_n*dbdp; 

    end

    Bt = -Bt;
    Bp = -1/sin(th)*Bp;

    % Compute sine and cosine of declination and right ascension
    cdec = cos(dec); sdec = sin(dec);
    cras = cos(ras); sras = sin(ras);

    % Compute B in terms of geocentric inertial components 
    Bx = (Br*cdec + Bt*sdec)*cras - Bp*sras; 
    By = (Br*cdec + Bt*sdec)*sras + Bp*cras; 
    Bz = Br*sdec - Bt*cdec;

    % Transform B from [nT] to [T]
    B = 1e-9*[Bx; By; Bz];
    
end

function [cth, sth] = sincos(maxdeg, th)

    % Efficiently compute the sin e cos of a given angle

    cth = zeros(maxdeg + 1, 1);
    sth = zeros(maxdeg + 1, 1);

    cth(1) = 1.0; 
    sth(1) = 0.0; 

    sm = sin(th); 
    cm = cos(th); 
    
    % Exploit the recursive sine and cosine relation to avoid repeated
    % calls to the trigonometric functions. 
    for m = 1:maxdeg
        cth(m+1) = cm*cth(m) - sm*sth(m);
        sth(m+1) = cm*sth(m) + sm*cth(m);
    end

end

function [P, dP] = compute_legendre(maxdeg, th)

    % Iteratively compute all the legendre polynomials and their derivatives 
    % in vector form up to maxdeg. 

    cth = cos(th); 
    sth = sin(th);

    % Vector size
    vsize = (maxdeg + 1)*(maxdeg + 2)/2;

    P = zeros(vsize, 1);
    dP = zeros(vsize, 1);

    P(1, 1) = 1.0;  
    dP(1, 1) = 0.0;

    for n = 1:maxdeg
        crtn = n*(n+1)/2 + 1;

        for m = 0:n
            idx = crtn + m;

            if m == n 
                prv = idx - n - 1; 
                P(idx) = sth*P(prv);
                dP(idx) = sth*dP(prv) + cth*P(prv);

            elseif m ~= n - 1
                prv = idx - n; 
                prvv = idx - 2*n + 1;

                Knm = ((n-1)^2 - m^2)/((2*n-1)*(2*n-3));

                P(idx) = cth*P(prv) - Knm*P(prvv);
                dP(idx) = cth*dP(prv) - sth*P(prv) - Knm*dP(prvv);

            else
                prv = idx - n; 
                P(idx) = cth*P(prv);
                dP(idx) = cth*dP(prv) - sth*P(prv);

            end
        end
    end

end
