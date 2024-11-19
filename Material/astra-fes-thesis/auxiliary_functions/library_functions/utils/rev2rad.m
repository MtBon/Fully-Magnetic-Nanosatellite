function rad = rev2rad(rev)

    % Convert revolutions to radians.
    %
    % Parameter
    % ---------
    %   rev: double 
    %         fraction of revolution, between [0, 1] 
    %   
    % Returns
    % -------
    %   rad: double 
    %         revolution angle in radians, between [0, 2π]
    
    rad = 2.0 * pi * (rev - fix(rev));

end
