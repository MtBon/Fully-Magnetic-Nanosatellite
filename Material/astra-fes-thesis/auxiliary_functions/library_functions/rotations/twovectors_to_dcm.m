
%{
    twovectors_to_dcm(a, b, seq)

Generate a Direction Cosine Matrix from two time-dependent vectors 
`a` and `b`, following the directions specified in `seq`. 

# Inputs 

- `a`: *number (3, N)* -- The primary vector that will be aligned with 
    the first direction specified in `seq`. 
- `b`: *number (3, N)* -- The secondary vector. The component of this 
    vector that is orthogonal to the primary vector is aligned with the 
    second direction specified in the sequence `seq`.
- `seq`: *string* -- Accepted sequence directions are: `XY`, `YX`, `XZ`,
    `ZX`, `YZ` or `ZY`.

# Outputs 
- `A`: *number (3, 3, N)* -- Desired Direction Cosine Matrix (DCM) that 
    rotates a vector defined in the initial references axes to those 
    formed by `a` and `b`.

# See also 
See also [`twovectors_to_ddcm`](@ref).

!!! note 
    The primary and secondary vectors do not have to be orthogonal.
    However, a great loss of precision happens when the two vectors 
    are almost aligned. This function does not perform any check on 
    the angular separation of the two vectors. The user should ensure 
    that the primary and secondary vector differ of at least 1 mrad.

%}

function A = twovectors_to_dcm(a, b, seq)

    %#codegen

    if seq == "XY" 
        w = cross(a, b);
        v = cross(w, a);
        u = a;
    
    elseif seq == "YX" 
        w = cross(b, a);
        u = cross(a, w);
        v = a;
    
    elseif seq == "XZ" 
        v = cross(b, a);
        w = cross(a, v);
        u = a ;
    
    elseif seq == "ZX" 
        v = cross(a, b);
        u = cross(v, a);
        w = a ;
    
    elseif seq == "YZ" 
        u = cross(a, b);
        w = cross(u, a);
        v = a;
    
    elseif seq == "ZY"
        u = cross(b, a);
        v = cross(a, u);
        w = a;
    else
        error("Invalid directions.")
    end

    % Compose all the rotation matrices
    A = zeros(3, 3);
    A(1, :) = vec3_to_unit(u); 
    A(2, :) = vec3_to_unit(v); 
    A(3, :) = vec3_to_unit(w);

end

function u = vec3_to_unit(x)

    u = x;
    n = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    if n ~= 0
        u = u/n;
    end

end
