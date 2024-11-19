function man = tan2man(theta, ecc)

    % Convert the true anomaly `tht` to the associated mean anomaly 
    % for an elliptical orbit of eccentricity `ecc`.

    %#codegen 
    
    if ecc == 0
        man = theta;
        
    elseif ecc < 1
        E = 2*atan(sqrt((1-ecc)./(1+ecc)).*tan(0.5*theta));
        man = E - ecc.*sin(E);
        
    elseif ecc == 1
        % Compute parabolic mean anomaly
        t = tan(0.5*theta);
        man = 0.5*t + t^3/6;

    else
        % Compute hyperbolic eccentric anomaly 
        F = 2*atanh(sqrt((ecc-1)/(ecc+1))*tan(0.5*theta));

        % Compute hyperbolic mean anomaly
        man = ecc*sinh(F) - F;
    end


end