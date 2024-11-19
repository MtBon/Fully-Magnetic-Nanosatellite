function R = dxddoe(coe, GM)

    % Compute the jacobian of the LVLH state with respect to the 
    % (follower - leader) orbital elements difference. 
    %
    % Parameters
    % ----------
    %   coe: double(6, 1)
    %       Classical Orbital Elements [sma, ecc, inc, ran, aop, tan]
    %   GM: double
    %       Planetary gravitational constant.
    %
    % Returns 
    % -------
    %   R: Jacobian of LVLH with respect to DOE.
    %
    % Notes
    % -----
    %   All the input quantities must have coherent units of measure
    %   (e.g., if GM is in km³/s², then coe(1) must be in km aswell).
    %
    % References 
    % ----------
    %   [1] D'Amico Simone, Autonomous Formation Flying in Low Earth Orbit,
    %       PhD Thesis, 2010
    %
    %   [2] Silvestrini Stefano, AI-Augmented Guidance, Navigation and 
    %       Control for Proximity Operations of Distributed Systems

    sma = coe(1); sma2 = sma*sma; sma3 = sma2*sma;
    ecc = coe(2);
    
    pi2 = 2*pi;
    
    % Wrap angles 
    inc = mod(coe(3), pi);
    aop = mod(coe(5), pi2);
    th  = mod(coe(6), pi2);

    stan = sin(th); ctan = cos(th);
    sinc = sin(inc); cinc = cos(inc); 
    saop = sin(aop); caop = cos(aop); 

    staop = stan*caop + ctan*saop;
    ctaop = ctan*caop - stan*saop;

    % Mean motion 
    n = sqrt(GM/sma3);
        
    % Preliminaries 
    eta = sqrt(1-ecc^2);
    r = sma*eta^2/(1+ecc*ctan); r2 = r*r;

    R11 = r/sma; 
    R12 = sma*ecc*stan/eta;
    R14 = -sma*ctan;

    R22 = sma2/r*eta;
    R23 = r; 
    R24 = (sma+r/eta^2)*stan;
    R26 = r*cinc;

    R35 = r*staop;
    R36 = -r*sinc*ctaop;
    
    R41 = -n*ecc*stan/(2*eta);
    R42 = ecc*n*ctan*sma3/r2;
    R44 = n*stan*eta*(sma3/r2);

    R52 = -ecc*n*stan*sma3/r2;
    R53 = sma*ecc*n*stan/eta;
    R54 = n*eta*(1+r/(sma*eta^2))*sma3/r2*ctan+sma*ecc*n*stan^2/eta^3;
    
    R65 = sma*n/eta*(ctaop+ecc*caop);
    R66 = sma*n/eta*(sinc*(staop+ecc*saop));
    
    R = [R11, R12,   0, R14,   0,   0;
           0, R22, R23, R24,   0, R26; 
           0,   0,   0,   0, R35, R36; 
         R41, R42,   0, R44,   0,   0; 
           0, R52, R53, R54,   0,   0; 
           0,   0,   0,   0, R65, R66];

end