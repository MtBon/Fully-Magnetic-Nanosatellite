function roe = coe2roe(kep_f, kep_r)

    % Compute relative orbital elements `roe` given the classical orbital 
    % elements of the follower `kep_f` and of the reference `kep_r`. 
    
    % The output vector equals: roe = [δa, δλ, δex, δey, δix, δiy]

    sma_f = kep_f(1); sma_r = kep_r(1); 
    ecc_f = kep_f(2); ecc_r = kep_r(2);
    inc_f = kep_f(3); inc_r = kep_r(3);
    ran_f = kep_f(4); ran_r = kep_r(4); 
    aop_f = kep_f(5); aop_r = kep_r(5); 
    tan_f = kep_f(6); tan_r = kep_r(6); 

    % Compute eccentricity vectors
    ex_f = ecc_f*cos(aop_f); 
    ex_r = ecc_r*cos(aop_r);

    ey_f = ecc_f*sin(aop_f); 
    ey_r = ecc_r*sin(aop_r); 

    % Compute mean anomalies 
    man_f = tan2man(tan_f, ecc_f); 
    man_r = tan2man(tan_r, ecc_r); 

    % Compute argument of latitudes 
    lat_f = aop_f + man_f; 
    lat_r = aop_r + man_r; 

    roe = [
            (sma_f - sma_r)/sma_r; 
            (lat_f - lat_r) + (ran_f - ran_r)*cos(inc_r);
            ex_f - ex_r; 
            ey_f - ey_r; 
            inc_f - inc_r; 
            (ran_f - ran_r)*sin(inc_r);
          ];

end


