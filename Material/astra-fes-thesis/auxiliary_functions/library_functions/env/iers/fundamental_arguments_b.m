function [Ma, Sa, um, Ds, Om] = fundamental_arguments_b(t)

    % Compute the Fundamental (Delaunay) Luni-Solar arguments, in radians, 
    % associated to the IAU2000B model at time `t`. 
    %
    % Parameters 
    % ----------
    %   t: double
    %       Time expressed in TDB Julian centuries since J2000
    %
    % Returns 
    % -------
    %   Ma: double
    %       Mean anomaly of the Moon [rad].
    %   Sa: double
    %       Mean anomaly of the Sun [rad].
    %   um: double
    %       Mean longitude of the Moon minus mean longitude of the ascending node [rad].
    %   Ds: double
    %       Mean elongation of the Moon from the Sun [rad].
    %   Om: double
    %       Mean longitude of the Moons ascending node [rad].
    %
    % References
    % ----------
    %   [1] Simon, J. et al., (1994), Numerical expressions for precession formulae
    %       and mean elements of the Moon and the planets.

    % Arcseconds in a full circle 
    ARCSECTURN = 1296000.0;

    % Arcseconds to radians conversion factor
    arcsec2rad = pi/648000;

	Ma = mod(485868.249036 + 1717915923.2178*t, ARCSECTURN) * arcsec2rad;
    
	Sa = mod(1287104.79305 + 129596581.0481*t, ARCSECTURN) * arcsec2rad;

	% Mean argument of the latitude of the Moon. 
	um = mod(335779.526232 + 1739527262.8478*t, ARCSECTURN) * arcsec2rad;

	% Mean elongation of the Moon from the Sun. 
	Ds = mod(1072260.70369 + 1602961601.2090*t, ARCSECTURN) * arcsec2rad;

	% Mean longitude of the ascending node of the Moon. 
	Om = mod(450160.398036 - 6962890.5431*t, ARCSECTURN) * arcsec2rad;

end