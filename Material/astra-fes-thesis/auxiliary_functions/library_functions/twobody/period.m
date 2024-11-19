function T = period(sma, GM)

    % Compute the period of an elliptical orbit given its 
    % semimajor axis `sma` and the gravitational parameter `GM`

    T = 2*pi*sqrt(sma^3/GM);
    
end
