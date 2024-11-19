function [dstv_lvlh] = car2lvlh(stv_ref, stv_fol, GM)

    % Obtain the relative (chaser-target) or (target-chaser) state in 
    % LVLH axes given the reference state vector `stv_t`, the follower 
    % state vector `stv_c` and the planetary constant `GM`.

    pos_ref = reshape(stv_ref(1:3), 3, 1); 
    vel_ref = reshape(stv_ref(4:6), 3, 1);

    pos_fol = reshape(stv_fol(1:3), 3, 1);
    vel_fol = reshape(stv_fol(4:6), 3, 1);

    % Compute target acceleration, assuming a keplerian motion
    % This assumption introduces an error in the velocity conversion. 
    acc_ref = -GM*pos_ref/norm(pos_ref)^3;

    % Compute rotation matrix from ECI to LVLH and its derivative
    A  = twovectors_to_dcm(pos_ref, vel_ref, "XY");
    dA = twovectors_to_ddcm([pos_ref; vel_ref], [vel_ref; acc_ref], "XY");
                    
    % Compute relative state in ECI 
    dpos_ECI = pos_fol - pos_ref;
    dvel_ECI = vel_fol - vel_ref;

    % Convert relative state from ECI to LVLH
    dpos_LVLH = A*dpos_ECI;
    dvel_LVLH = A*dvel_ECI + dA*dpos_ECI;

    dstv_lvlh = [dpos_LVLH; dvel_LVLH];

end