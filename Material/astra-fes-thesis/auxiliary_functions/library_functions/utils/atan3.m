function y = atan3(a, b)

    % Four quadrant inverse tangent of a and b. 
    % 
    % Parameters
    % ----------
    %   a: double 
    %       sine of the angle
    %   b: double 
    %       cosine of the angle
    %
    % Returns
    % -------
    %   y: double 
    %       output angle (rad) between [0, 2π]

    epsilon = 0.0000000001;
    pidiv2 = 0.5 * pi;

    if (abs(a) < epsilon)
        
        if (abs(b) < epsilon)
            y = 0; 
        else
            y = (1 - sign(b)) * pidiv2;
        end
        
        return;
    end

    c = (2 - sign(a)) * pidiv2;

    if (abs(b) < epsilon)
        y = c;
    else
        y = c + sign(a) * sign(b) * (abs(atan(a / b)) - pidiv2);
    end

end