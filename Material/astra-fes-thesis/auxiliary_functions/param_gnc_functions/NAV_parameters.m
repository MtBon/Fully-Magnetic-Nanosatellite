function [TUN_STRUCT, NTUN_STRUCT, SAMPLE_TIME] = NAV_parameters

    % It returns a struct array containing tunable struct parameters
    % for NAV Block Models
    
    % Retrieve the spacecraft initial state in ECI (in m and m/s)
    % X0_sc = evalin('caller', 'x_sc');

    % Retrieve Gravitational parameter (m³/s²)
    GM_EARTH = evalin('caller', 'env.GM_EARTH');
    
    % Sample time (s)
    SAMPLE_TIME.time_1 = 1;
    % ... add other sample times if present
    
    % --- TUNABLE Parameters

    % IB to ICB offsets and DCM
    TUN_STRUCT.offset_IB2ICB = single(zeros(3, 1));  % TODO: FIXAMI 
    TUN_STRUCT.A_ICB2IB = single(eye(3));            % TODO: FIXAMI

    % --- NON TUNABLE Parameters
    % Earth planetary constant (m³/s²)
    NTUN_STRUCT.GM_EARTH = single(GM_EARTH); 

    % Earth radius (m)
    NTUN_STRUCT.RADIUS_EARTH = single(evalin('caller', 'env.RADIUS_EARTH'));
    
    % Earth J2 Coefficient (-)
    NTUN_STRUCT.J2_EARTH = single(evalin('caller', 'env.J2'));


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If Needed GNC parameters for sensors and actuators are contained in their
    % init functions
    %
    % Example to Initialize Actuators GNC Data  
    % ACT struct for data used in GNC - To be used in eventual Param GNC Function
    % ACT_DKE struct for data used in DKE
    % [MTQ,~]=HP_MT_GST600_init(MTQ_actuator_time); %#ok<ASGLU> - are needed in the model of MTQ - if MTQ are present, MTQ_DKE shall be present 
    % [RWL,~]=HP_RW_GSW600_init(RWL_actuator_time); %#ok<ASGLU> - are needed in the dynamics and in the model of RWL - RWL_DKE shall be always present

end
