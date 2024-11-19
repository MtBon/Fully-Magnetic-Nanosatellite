function R = angle1_to_dcm(x, rot_axis)

    % Create a direction cosine matrix that performs a rotation x about the coordinate 
    % axis specified in `rot_axis`
    %
    % Parameters 
    % ----------
    %   x: double 
    %       Rotation angle 
    %   rot_axis: string
    %       Rotation axis. The possible values are: 'X', 'Y' or 'Z'
    %  
    % Returns 
    % -------
    %   R: double(3, 3)
    %       Desired rotation matrix.

    c = cos(x); s = sin(x); 

    R = eye(3);

    if strcmp(rot_axis, "X")
        R = [1,  0, 0; 
             0,  c, s; 
             0, -s, c];

    elseif strcmp(rot_axis, "Y")
        R = [c,  0, -s; 
             0,  1,  0; 
             s,  0,  c]; 
         
    elseif strcmp(rot_axis, "Z")
        R = [ c, s, 0; 
             -s, c, 0; 
              0, 0, 1];
    end

end