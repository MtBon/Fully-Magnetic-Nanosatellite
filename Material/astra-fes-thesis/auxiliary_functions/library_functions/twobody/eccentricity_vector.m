
function e = eccentricity_vector(stv, GM)

    % Compute the eccentricity vector given the position `r`, velocity `v` and the 
    % planetary constant `GM`.
    
    r = stv(1:3); 
    v = stv(4:6);

    e = ((dot(v, v) - GM/norm(r))*r - v*dot(r, v))/GM;
    
end
