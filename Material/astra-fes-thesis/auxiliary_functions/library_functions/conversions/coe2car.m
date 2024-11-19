
%{
    coe2car(coe, GM)

Convert the Classical Orbital Elements (COE) into an 
inertial cartesian state (position and velocity).

# Inputs 
- `coe`: *number (6, N)* -- Classical orbital elements. The elements
    of each column vector are organised as follows: Semi-major axis, 
    eccentricity, inclination, right ascension of the ascending node, 
    argument of perigee, true anomaly. 
- `GM`: *number* -- Planetary gravitational constant. 

# Outputs 
- `stv`: *number (6, N)* -- Cartesian inertial state vector (position and velocity).

# See Also 
See also [`kat_car2coe`](@ref).

# References 
- D. Vallado, _Fundamentals of Astrodynamics_, IV Edition, 2013

!!! note 
    All the angular quantities within `coe` must be expressed in radians. 
    The semi-major axis must be expressed in the same units of length of `GM`.

%}

function stv = coe2car(coe, GM)

    %#codegen
    sma = coe(1); ecc = coe(2); inc = coe(3); 
    ran = coe(4); aop = coe(5); theta = coe(6); 

    ctan = cos(theta); 

    fpa = atan2(ecc.*sin(theta), 1 + ecc.*ctan);

    lan = aop + theta; 
    flan = lan - fpa;

    r = sma.*(1 - ecc.^2)./(1 + ecc.*ctan);
    v = sqrt(2*GM./r - GM./sma);

    sinc = sin(inc);  cinc = cos(inc);  
    sran = sin(ran);  cran = cos(ran); 
    slan = sin(lan);  clan = cos(lan);
    sfla = sin(flan); cfla = cos(flan);

    cs1 = slan.*cinc;
    cs2 = cfla.*cinc;

    stv = [r.*(clan.*cran - cs1.*sran); 
           r.*(clan.*sran + cs1.*cran); 
           r.*slan.*sinc;
           v.*(-sfla.*cran - cs2.*sran);
           v.*(-sfla.*sran + cs2.*cran); 
           v.*cfla.*sinc];
       
end