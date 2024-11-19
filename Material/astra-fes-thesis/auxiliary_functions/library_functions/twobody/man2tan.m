
function theta = man2tan(man, ecc)
    
    %#codegen
    if ecc == 0
        theta = man; 
        
    elseif ecc < 1
        theta = man2tan_elliptical(man, ecc); 
    
    elseif ecc == 1
        theta = man2tan_parabolic(man); 

    else
        theta = man2tan_hyperbolic(man, ecc);

    end

end

function theta = man2tan_elliptical(man, ecc)

    tol = 1e-9;
    err = 1 + tol;
    if man < pi 
        Ep = man + 0.5*ecc;
    else
        Ep = man - 0.5*ecc;
    end

    i = 0;
    while (err > tol) && (i < 100)

        f = ecc*sin(Ep) - Ep + man;
        df = ecc*cos(Ep) - 1;

        E = Ep - f/df;

        err = abs(E - Ep);
        Ep = E;
        i = i + 1;
    end

    theta = mod(2*atan(sqrt((1+ecc)/(1-ecc))*tan(0.5*E)), 2*pi);

end

function theta = man2tan_parabolic(man)

    z = (3*man + sqrt(1 + 9*man.^2)).^(1/3);
    theta = 2*atan(z - 1/z);

end

function theta = man2tan_hyperbolic(man, ecc)

    tol = 1e-9;
    err = 1 + tol;

    % Initial guess is roughly assumed equal.
    Fp = man;

    i = 0;
    while (err > tol) && (i < 100)

        f = ecc*sinh(Fp) - Fp - man;
        df = ecc*cosh(Fp) - 1;
        
        F = Fp - f/df; 

        err = abs(F - Fp); 
        Fp = F; 
        i = i + 1;
    end

    theta = 2*atan(sqrt((ecc+1)/(ecc-1))*tanh(0.5*F));

end