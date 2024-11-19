function C_x2roe = hill2roe_1storder(coe_r, GM)

    % From relative orbital elements to hill frame of the reference satellite
    % using 1st order approximation. 
    
    % Conversion is obtained as: dROE/dX =  dROE/dDOE * dDOE/dX
    
    sma = coe_r(1); ecc = coe_r(2);
    inc = coe_r(3); aop = coe_r(5);

    sinc = sin(inc); cinc = cos(inc); 
    saop = sin(aop); caop = cos(aop); 
    
    % dROE/dDOE
    dROEdDOE = [1/sma,  0,          0,     0,  0,    0;
                    0,  1,          1,     0,  0, cinc;
                    0,  0,  -ecc*saop,  caop,  0,    0; 
                    0,  0,   ecc*caop,  saop,  0,    0; 
                    0,  0,          0,     0,  1,    0; 
                    0,  0,          0,     0,  0, sinc];
    
    % Compute matrix
    dXdDOE = dxddoe(coe_r, GM);
    C_x2roe = (dROEdDOE / dXdDOE);
    
end
