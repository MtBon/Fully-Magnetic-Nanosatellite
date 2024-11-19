function roe = get_starting_point(idx, sma_t, n)

    % Return the relative orbital elements `roe` of a given holding point/
    % orbit or inspection orbit.
    % Inputs:
    % - number `idx`;
    % - target semi-major axis `sma_t` [m];
    % - target mean motion 'n';
    
    switch idx
        case 1
            % Hold point
            droe = [0; 20e3; 0; 0; 0; 0];

        case 2 
            % Hold orbit
            droe = [0 -5e3 0 1e3 0 1e3];
            
        case 3 
            % Hold orbit
            droe = [0 2e3 0 500 0 500];

        case 4
            % Hold orbit 
            droe = [0 -2e3 0 200 0 200];

        case 5
            % Drifting ...
            droe = [25000/(24*3600*1.5*n) 20e3 0 1e3 0 1e3];

        case 6 
            % Drifting ...
            droe = [-7000/(0.5*24*3600*1.5*n) -5e3 0 500 0 500];
            
        case 7 
            droe = [4000/(0.5*24*3600*1.5*n) 2e3 0 200 0 200];

        case 8
            droe = [-4000/(1.5*24*3600*1.5*n) -2e3 0 100 0 100];

        otherwise 
            error("Unrecognised Holding Point");
    end

    roe = droe/sma_t; 

end