function R = angle3_to_dcm(x, y, z, rot_seq)

    % Create a direction cosine matrix that performs a set of 
    % rotations (x, y, z) about the coordinate axes specified in `rot_seq`
    %
    % Parameters 
    % ----------
    %   x: double 
    %       First rotation angle 
    %   y: double
    %       Second rotation angle
    %   z: double
    %       Third rotation angle
    %   rot_seq: string
    %       Rotation sequence. The possible values are: `XYX`, `XYZ`, `XZX`, `XZY`, 
    %       `YXY`, `YXZ`, `YZX`, `YZY`, `ZXY`, `ZXZ`, `ZYX`, or `ZYZ`
    %  
    % Returns 
    % -------
    %   R: double(3, 3)
    %       Desired rotation matrix.

    c1 = cos(x); s1 = sin(x); 
    c2 = cos(y); s2 = sin(y); 
    c3 = cos(z); s3 = sin(z);

    R = eye(3);

    if strcmp(rot_seq, "ZYX")
        R = [           c2*c1,            c2*s1,    -s2;
             s3*s2*c1 - c3*s1, s3*s2*s1 + c3*c1,  s3*c2;
             c3*s2*c1 + s3*s1, c3*s2*s1 - s3*c1,  c3*c2];

    elseif strcmp(rot_seq, "XYX")
        R = [   c2,             s1*s2,           -c1*s2;
             s2*s3, -s1*c2*s3 + c1*c3, c1*c2*s3 + s1*c3;
             s2*c3, -s1*c3*c2 - c1*s3, c1*c3*c2 - s1*s3];

    elseif strcmp(rot_seq, "XYZ")
        R = [c2*c3,  s1*s2*c3 + c1*s3, -c1*s2*c3 + s1*s3;
            -c2*s3, -s1*s2*s3 + c1*c3,  c1*s2*s3 + s1*c3;
                s2,            -s1*c2,             c1*c2];

    elseif strcmp(rot_seq, "XZX")
        R = [   c2,             c1*s2,             s1*s2;
            -s2*c3,  c1*c3*c2 - s1*s3,  s1*c3*c2 + c1*s3;
             s2*s3, -c1*c2*s3 - s1*c3, -s1*c2*s3 + c1*c3];

    elseif strcmp(rot_seq, "XZY")
        R = [c3*c2, c1*c3*s2 + s1*s3, s1*c3*s2 - c1*s3;
               -s2,            c1*c2,            s1*c2;
             s3*c2, c1*s2*s3 - s1*c3, s1*s2*s3 + c1*c3];

    elseif strcmp(rot_seq, "YXY")
        R = [-s1*c2*s3 + c1*c3,  s2*s3, -c1*c2*s3 - s1*c3;
                         s1*s2,     c2,             c1*s2;
              s1*c3*c2 + c1*s3, -s2*c3,  c1*c3*c2 - s1*s3];

    elseif strcmp(rot_seq, "YXZ")
        R = [c1*c3 + s2*s1*s3, c2*s3, -s1*c3 + s2*c1*s3;
            -c1*s3 + s2*s1*c3, c2*c3,  s1*s3 + s2*c1*c3;
                        s1*c2,   -s2,             c2*c1];

    elseif strcmp(rot_seq, "YZX")
        R = [           c1*c2,     s2,            -s1*c2;
            -c3*c1*s2 + s3*s1,  c2*c3,  c3*s1*s2 + s3*c1;
             s3*c1*s2 + c3*s1, -s3*c2, -s3*s1*s2 + c3*c1];

    elseif strcmp(rot_seq, "YZY")
        R = [c1*c3*c2 - s1*s3, s2*c3, -s1*c3*c2 - c1*s3;
                       -c1*s2,    c2,             s1*s2;
             c1*c2*s3 + s1*c3, s2*s3, -s1*c2*s3 + c1*c3];

    elseif strcmp(rot_seq, "ZXY")
        R = [c3*c1 - s2*s3*s1, c3*s1 + s2*s3*c1, -s3*c2;
                       -c2*s1,            c2*c1,     s2;
             s3*c1 + s2*c3*s1, s3*s1 - s2*c3*c1,  c2*c3];

    elseif strcmp(rot_seq, "ZXZ")
        R = [-s1*c2*s3 + c1*c3, c1*c2*s3 + s1*c3, s2*s3;
             -s1*c3*c2 - c1*s3, c1*c3*c2 - s1*s3, s2*c3;
                         s1*s2,           -c1*s2,    c2];

    elseif strcmp(rot_seq, "ZYZ")
        R = [c1*c3*c2 - s1*s3,  s1*c3*c2 + c1*s3, -s2*c3;
            -c1*c2*s3 - s1*c3, -s1*c2*s3 + c1*c3,  s2*s3;
                        c1*s2,             s1*s2,     c2];
    end

end