function Q = xys2m(x, y, s)
    
    % Compute the Intermediate-to-Celestial matrix Q
    %
    % Parameters
    % ----------
    %   x: double 
    %       CIP x-coordinate [rad]. 
    %   y: double 
    %       CIP y-coordinate [rad]. 
    %   s: double
    %       CIO Locator [rad]
    % 
    % Returns
    % -------
    %   Q: double(3, 3)
    %       Intermediate-to-Celestial rotation matrix
    %
    % References
    % ----------
    %   [1] Wallace P. T. and Capitaine N. (2006), Precession-nutation
    %   procedures consistent with the IAU 2006 resolutions 
    
    % Retrieve the spherical angles E and d    
    r2 = x^2 + y^2; 
    E = atan3(y, x); 
    d = atan3(sqrt(r2), sqrt(1.0 - r2));
    
    Q = angle3_to_dcm(E+s, -d, -E, "ZYZ");

end