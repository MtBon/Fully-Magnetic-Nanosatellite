

%{
    mean_motion(sma, GM)

Compute the mean angular motion of a 2-Body orbit given its semi-major
axis `sma` and the planetary gravitational constant `GM`.

# See Also 
See also [`orbital_period`](@ref).

!!! note 
    It is assumed that `sma` and `GM` are expressed in the same units
    of measure of length. The returned value is measured in the same 
    units of time in which `GM` is expressed.

!!! warning 
    This function is valid only for circular and elliptical 
    orbits. 
%}

function T = mean_motion(sma, GM)

    %#codegen
    T = sqrt(GM./sma.^3);

end