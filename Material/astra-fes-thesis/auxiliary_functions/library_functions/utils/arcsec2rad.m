function y = arcsec2rad(x)
    
    % Convert arcseconds to radians. 
    %
    % Parameter
    % ---------
    %   x: double 
    %       Radians
    %  
    % Returns
    % -------
    %   y: double 
    %       Arcseconds
    
    y = x*pi/648000;

end