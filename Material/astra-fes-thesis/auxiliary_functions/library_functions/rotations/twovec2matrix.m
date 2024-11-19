function M=twovec2matrix(VEC,TARG)
%% Function TwoVec2Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Find rotation matrix from a given vector VEC to a given vector TARG
%
% © MAT - Aerospace Science and Technology Dept. - PoliMi - 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the rotation, vectors must be normalized
M=zeros(3,3);

VECn = VEC./norm(VEC);
TARGn = TARG./norm(TARG);

if any([any(isnan(VECn)),any(isnan(TARGn))])
    errax=[0;0;1]; %#ok<NASGU>
    errang=0; %#ok<NASGU>
else
    errax = cross(VECn, TARGn);
    errax = errax./norm(errax);
    % if cross(VECn, TARGn) is zero, vectors are parallel (angle = 0) or antiparallel
    % (angle = pi). In both cases it is necessary to provide a valid axis. Let's
    % select one that satisfies both cases - an axis that is perpendicular to
    % both vectors. We find this vector by cross product of the first vector
    % with the "least aligned" basis vector.
    if ~any(errax)
        absVEC = abs(VECn);
        [~, mind] = min(absVEC);
        c = zeros(3,1);
        c(mind) = 1;
        errax = cross(VECn, c);
        errax=errax./norm(errax);
    end
    % min to eliminate possible rounding errors that can lead to dot product >|1|
    errang = acos( max(min(dot(VECn, TARGn), 1), -1) );
    
    % build the rotation matrix
    s = sin(errang);
    c = cos(errang);
    t = 1 - c;
    
    x = errax(1);
    y = errax(2);
    z = errax(3);
    M = [ ...
        t*x*x + c,    t*x*y - s*z,  t*x*z + s*y; ...
        t*x*y + s*z,  t*y*y + c,    t*y*z - s*x; ...
        t*x*z - s*y,  t*y*z + s*x,  t*z*z + c ...
        ];
end
end