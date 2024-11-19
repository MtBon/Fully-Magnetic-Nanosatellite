function acc = gravity_field_spherical(pos, GM, R, Shn, Chn, Khn, maxdeg)

    % Compute the perturbing gravitational acceleration in the body-fixed
    % axes using the gravitational spherical harmonics expansion. The 
    % output acceleration does NOT include the zero-order spherical
    % contribution.
    %
    % Parameters
    % ----------
    %   pos: double (3, 1)
    %       Cartesian position vector in the body-fixed frame in [m]. 
    %   GM: double
    %       Gravitational parameter of the gravity model in [m³/s²]
    %   R:  double 
    %       Reference body radius of the gravity model in [m]
    %   Shn: double (N, 1)
    %       Vector storing the normalised sine coefficients computed from
    %       the parse_grav_data.m function. 
    %   Chn: double (N, 1)
    %       Vector storing the normalised cosine coefficients computed from
    %       the parse_grav_data.m function. 
    %   Khn: double (N, 1)
    %       Legendre polynomial scale factors.
    %   max_deg: integer 
    %       Max degree and order of the harmonic expansion.
    %
    % Returns 
    % -------
    %   acc: double (3, 1)
    %       Gravitational acceleration in the body-fixed frame in [m/s²].
    %
    % Notes
    % -----
    %   The ScaleFct of the Legendre's polynomials is requested to properly 
    %   dimensionalise the dudlat equation (see [1] Efficient models for 
    %   the evaluation and estimation of the gravity field, B. Jones, PhD Thesis). 
    %
    %   The classic algorithm has a singularity at the poles because the
    %   derivative of drdp has (x^2 + y^2) at the denominator. Possible
    %   workarounds are the Pines or Cartesian formulation (however [2] says
    %   the Pines formulation suffers from intrisic loss of precision at the
    %   equator. For now, this algorithm implements the MATLAB solution.
    %
    %   See ([2] Fantino E., Casotto S., Methods of Harmonic Synthesis for 
    %   Global Geopotential Models and Their First-, Second-, and Third-order
    %   Gradients) for a detailed description of a new efficient way to
    %   implement this algorithm, including comparisons with other models.
    
    x = pos(1); y = pos(2); z = pos(3);

    % Convert Body-fixed position to bodycentric latitude and longitude.
    [lat, lon, ~] = pos2geoc(pos);

    xy2 = x^2 + y^2; 
    r2 = xy2 + z^2; 
    r = sqrt(r2);

    % Evaluates normalised Legendre's polynomials
    P = compute_legendre(maxdeg, lat);

    % Pre-compute sine and cosine of longitude iteratively
    [clon, slon] = sincos(maxdeg, lon);

    % Precompute tan of latitude 
    tlat = tan(lat);

    r_ratio = R/r;
    r_ratio_j = r_ratio;
    
    dudr = 0.0; % Excludes the spherical gravitational contribution!
    dudp = 0.0;
    dudl = 0.0; 

    for j = 2:maxdeg
        IDx = 0.5*j*(j+1) + 1;
        
        dudr_k = 0.0; 
        dudp_k = 0.0; 
        dudl_k = 0.0;
        
        for k = 0:j
            
            % Extract coefficients 
            Cjk = Chn(IDx+k);
            Sjk = Shn(IDx+k);
            Pjk = P(IDx+k);

            slk = slon(k+1);
            clk = clon(k+1);
            
            csjk = Cjk*clk + Sjk*slk;

            % Radial partial derivative
            dudr_k = dudr_k + Pjk*csjk;
            
            % Latitude partial derivative
            if k == j 
                % In this case nPjk = 0.0! 
                dudp_k = dudp_k - k*tlat*Pjk*csjk;
            else 
                Kjk = Khn(IDx+k);
                nPjk = P(IDx+k+1);
                dudp_k = dudp_k + (nPjk*Kjk - k*tlat*Pjk)*csjk;
            end

            % Longitude partial derivative
            if k ~= 0
                dudl_k = dudl_k + k*Pjk*(Sjk*clk - Cjk*slk);
            end
        end
        
        r_ratio_j = r_ratio_j*r_ratio;

        dudr = dudr + (j+1)*r_ratio_j*dudr_k;
        dudp = dudp + r_ratio_j*dudp_k; 
        dudl = dudl + r_ratio_j*dudl_k; 

    end

    tmp = GM/r; 

    dudr = -tmp/r*dudr; 
    dudp = tmp*dudp; 
    dudl = tmp*dudl;

    xy = sqrt(xy2); 

    % Handle the singularity at the poles 
    if pi/2 - abs(lat) <= 1e-8
        gx = 0; 
        gy = 0; 
        gz = dudr*z/r;
        
    else 
        e1 = dudr/r;
        e2 = (e1 - z/(r2*xy)*dudp);
        e3 = dudl/xy2; 

        % Transform acceleration from spherical to cartesian body-fixed
        % coordinates.
        gx = e2*x - e3*y;
        gy = e2*y + e3*x;
        gz = e1*z + xy/r2*dudp;
    end

    acc = [gx; gy; gz];

end

function [clon, slon] = sincos(maxdeg, lon)

    % Efficiently compute the sin e cos of a given angle
    clon = zeros(maxdeg+2, 1);
    slon = zeros(maxdeg+2, 1);

    clon(1) = 1.0; 
    slon(1) = 0.0; 

    sl = sin(lon); 
    cl = cos(lon); 
    
    % Exploit the recursive sine and cosine relation to avoid repeated
    % calls to the trigonometric functions. 
    for m = 1:maxdeg+1
        clon(m+1) = cl*clon(m) - sl*slon(m);
        slon(m+1) = cl*slon(m) + sl*clon(m);
    end

end

function P = compute_legendre(maxdeg, lat)

    % P and S must be vectors of (max_deg+3)*(max_deg+4)/2
    vsize = (maxdeg+3)*(maxdeg+4)/2;

    P = zeros(vsize, 1);

    x = sin(lat); 
    y = cos(lat);

    P(1) = 1.;
    P(2) = sqrt(3)*x;
    P(3) = sqrt(3)*y;

    for j = 2:maxdeg+2

        IDx = j + 1;

        crt = IDx*(IDx+1)/2;
        prv = j*IDx/2;
        pprv = j*(IDx-2)/2;

        P(crt) = sqrt(1.0+0.5/j)*y*P(prv);

        for k = 0:j-1
            if k == 0
                P(crt-j+k) = sqrt(2*j+1)/j*(sqrt(2*j-1)*x*P(prv-j+1+k) - (j-1)/sqrt(2*j-3)*P(pprv-j+2+k));
            else 
                P(crt-j+k) = sqrt(2*j+1)/sqrt((j+k)*(j-k))*(sqrt(2*j-1)*x*P(prv-j+1+k) - sqrt(j+k-1)*sqrt(j-k-1)/sqrt(2*j-3)*P(pprv-j+2+k));
            end
        end

    end

end