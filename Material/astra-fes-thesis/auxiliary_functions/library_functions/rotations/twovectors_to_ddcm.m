
%{
    twovectors_to_ddcm(a, b, seq)

Generate the time derivative of a Direction Cosine Matrix from two
time-dependent state vectors `a` and `b`, following the directions 
specified in `seq`. 

# Inputs 

- `a`: *number (6, N)* -- The primary vector that will be aligned with 
    the first direction specified in `seq`. 
- `b`: *number (6, N)* -- The secondary vector. The component of this 
    vector that is orthogonal to the primary vector is aligned with the 
    second direction specified in the sequence `seq`.
- `seq`: *string* -- Accepted sequence directions are: `XY`, `YX`, `XZ`,
    `ZX`, `YZ` or `ZY`.

# Outputs 
- `A`: *number (3, 3, N)* -- Desired Direction Cosine Matrix (DCM) that 
    rotates a vector defined in the initial references axes to those 
    formed by `a` and `b`.

# See also 
See also [`twovectors_to_dcm`](@ref).

!!! note 
    The input state vectors must have in the 4th to 6th rows the 
    derivatives of the elements in the first three rows, respectively.

%}

function A = twovectors_to_ddcm(a, b, seq)
    
    %#codegen

    if seq == "XY" 
        w = cross6(a, b);
        v = cross6(w, a);
        u = a;
    
    elseif seq == "YX" 
        w = cross6(b, a);
        u = cross6(a, w);
        v = a;
    
    elseif seq == "XZ" 
        v = cross6(b, a);
        w = cross6(a, v);
        u = a ;
    
    elseif seq == "ZX" 
        v = cross6(a, b);
        u = cross6(v, a);
        w = a ;
    
    elseif seq == "YZ" 
        u = cross6(a, b);
        w = cross6(u, a);
        v = a;
    
    elseif seq == "ZY"
        u = cross6(b, a);
        v = cross6(a, u);
        w = a;
    else
        error("Invalid directions.")
    end

    % Compose all the rotation matrices
    A = zeros(3, 3);
    A(1, :) = vec3_to_dunit(u); 
    A(2, :) = vec3_to_dunit(v); 
    A(3, :) = vec3_to_dunit(w);

end

function w = cross6(u, v)

    % Compute the cross product of two vectors and its time derivative.
    %
    % Parameters
    % ----------
    %   u: number(6, N)
    %       First vector. If a matrix is provided, each column is treated 
    %       as a separate vector.
    %   v: number(6, N)
    %       Second vector. If a matrix is provided, each column is treated 
    %       as a separate vector.
    % 
    % Returns 
    % -------
    %   w: number(6, N)
    %       Output vector.
    %
    % Notes
    % -----
    %   The last elements of the inputs vectors must represent the 
    %   time derivatives of the first three.

    %#codegen 
    
    w = [u(2).*v(3) - u(3).*v(2); 
         u(3).*v(1) - u(1).*v(3); 
         u(1).*v(2) - u(2).*v(1); 
         u(5).*v(3) + u(2).*v(6) - u(6).*v(2) - u(3).*v(5);
         u(6).*v(1) + u(3).*v(4) - u(4).*v(3) - u(1).*v(6);
         u(4).*v(2) + u(1).*v(5) - u(5).*v(1) - u(2).*v(4)];

end

function y = vec3_to_dunit(x)
    
    % Compute the time derivative of a unit vector x.
    %
    % Parameters
    % ----------
    %   x: number(6, N)
    %       Input 6-elements vector. The elements in the 4th to 6th 
    %       positions must be the time derivatives of the respective
    %       elements in the 1st to 3rd positions. If a matrix is provided, 
    %       each column is treated as a separate vector.
    %   
    % Returns 
    % -------
    %   y: double(3, N)
    %       Unit vector derivative.
    %
    % Notes
    % -----
    %   The input vector can have arbitrary length, however only the 
    %   first 6 elements will be considered to compute the norm.

    %#codegen
    
    % Compute the vector norm
    r2 = x(1).*x(1) + x(2).*x(2) + x(3).*x(3);
    r = sqrt(r2); 
    r3 = r2.*r;

    d = -(x(1).*x(4) + x(2).*x(5) + x(3).*x(6))./r3;
    
    y = [x(4)./r + d.*x(1); 
         x(5)./r + d.*x(2); 
         x(6)./r + d.*x(3)];

end
