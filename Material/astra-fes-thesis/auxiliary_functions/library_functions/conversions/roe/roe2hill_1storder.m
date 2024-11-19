function C_roe2x = roe2hill_1storder(coe_r, GM)
    
    % Compute the transformation matrix from Relative Orbital Elements 
    % (ROE) to the cartesian LVLH state, using a 1st order approximation of
    % the transformation. 
    %
    % Parameters
    % ----------
    %   coe_r: double(6, 1)
    %       Classical Orbital Elements of the reference spacecraft, 
    %       [sma, ecc, inc, ran, aop, tan]
    %   GM: double
    %       Planetary gravitational constant 
    %
    % Returns 
    % -------
    %   C_roe2x: double(6, 6)
    %       Transformation matrix.
    %
    % References
    % ----------
    % 
    
    % Conversion is obtained as: dX/dROE =  dDOE/dROE * dX/dDOE

    sma = coe_r(1); ecc = coe_r(2);
    inc = coe_r(3); aop = coe_r(5);

    sinc = sin(inc); cinc = cos(inc); 
    saop = sin(aop); caop = cos(aop); 

    % dDOE/dROE 
    dDOEdROE = [sma,  0,         0,         0,  0,          0;
                  0,  1,  saop/ecc, -caop/ecc,  0, -cinc/sinc
                  0,  0, -saop/ecc,  caop/ecc,  0,          0;
                  0,  0,      caop,      saop,  0,          0;
                  0,  0,         0,         0,  1,          0; 
                  0,  0,         0,         0,  0,     1/sinc];

    % Compute matrix
    dXdDOE = dxddoe(coe_r, GM);
    C_roe2x = (dXdDOE*dDOEdROE);
    
end
