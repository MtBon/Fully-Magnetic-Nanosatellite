function coe_f = roe2coe(roe, coe_r)

    % Compute the classical orbital elements of follower `coe_f` 
    % from the relative orbital elements `roe` and the keplerian 
    % elements of the reference `coe_t`.

    % The relative orbital elements are written as: roe = [δa, δλ, δex, δey, δix, δiy]

    pi2 = 2*pi; 
    
    sma_r = coe_r(1); ecc_r = coe_r(2); inc_r = coe_r(3);
    ran_r = coe_r(4); aop_r = coe_r(5); tan_r = coe_r(6); 

    sma_f = sma_r*(roe(1) + 1); 

    % Retrieve eccentricity vectors
    ex_r = ecc_r*cos(aop_r);
    ey_r = ecc_r*sin(aop_r);

    ex_f = ex_r + roe(3); 
    ey_f = ey_r + roe(4);

    aop_f = mod(atan2(ey_f, ex_f), pi2); 
    saop_f = sin(aop_f); 
    caop_f = cos(aop_f);

    % Avoid division by 0.0
    if saop_f >= caop_f
        ecc_f = ey_f/saop_f; 
    else
        ecc_f = ex_f/caop_f;
    end

    cinc_r = cos(inc_r); 
    sinc_r = sin(inc_r);
        
    % Compute inclination and right ascension
    inc_f = mod(roe(5) + inc_r, pi);
    ran_f = mod(ran_r + roe(6)/sinc_r, pi2);

    % Retrieve chaser true anomaly.
    man_r = tan2man(tan_r, ecc_r);
    lan_r = aop_r + man_r;
    lan_f = roe(2) - cinc_r*(ran_f - ran_r) + lan_r;

    man_f = lan_f - aop_f;
    tan_f = man2tan(man_f, ecc_f);

    coe_f = [sma_f; ecc_f; inc_f; ran_f; aop_f; tan_f];

end