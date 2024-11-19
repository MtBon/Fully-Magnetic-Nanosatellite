close all
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 6DOF Dynamics - Orbit Attitude with GNC and Environment Simulator -
% Initialization file
%
% Full GNC/ADCS with AODS - ACS - Sensors - Actuators and Post Processing
%
%
% (c) ASTRA - Aerospace Science and Technology Dept. - PoliMi - 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Changelog:
%
%  V1.0 - HERMES ADCS FES - Andrea Colagrossi - < 2023
%  V2.0 - ASTRA FES - Andrea Colagrossi/Stefano Silvestrini - 16/01/2023
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../auxiliary_functions'))
addpath(genpath('../input_files'))
addpath(genpath('../assets'))

warning('off', 'Simulink:Engine:UINotUpdatedDuringRapidAccelSim');

%% SIMULATION TYPES AND INPUTS

% Requires the gnc_mode_flag to be simulated:
% 01 - Detumbling - DET
% 02 - Sun Pointing - SUN
% xx - Specific System Mode

% -- SIMULATION TYPE
% Number of simulation runs
n_sim = 1;

% True for repeatable Monte Carlo simulations
repeatable = 0;

% -- SIMULATION OPTIONS
% Defines the type of DKE I.C. simulation (e.g., detumbling, pointing, far
% range, ....)
sim_mode_type = 'pointing'; 

% Includes flexibility in the attitude dynamics
flex_flag = 0; 

% If TRUE, includes determination/navigation and sensors performance model
% disturbances (ONLY for the performance models)
gnc_pm_disturbance_flag = 1; 

% If TRUE, includes actuators performance model disturbances
acts_disturbance_flag = 1; 

% -- GNC MODE
% The value goes accordingly to the mode management definition
gnc_mode_flag = 1;

% -- OUTPUT FILE NAME
out_name = 'SimOutput';

% Select SIMULINK Model
if flex_flag
    % The flexible model is still not available in the DKE library
    error('FlexibleSimulator - v2.0 - Not Available Yet');
else
    simulator = 'ASTRA_FES_Simulator.slx';
end

%% RANDOM NUMBER SETTINGS

if repeatable
    % Rememeber to change seed value if want to repeat different simulations
    seed = 1;
    rng(seed);
else
    rng('shuffle')
end

%% FAILURE SIMULATOR 

% Variance of multiplier (Mean = 1) to simulate high noise of broken actuator
Fail_variance = 0.4/3;

%% ENVIRONMENT Input Files
% Read Input Files: Scenario definition - Reference Orbit - TARGET
filename_time = 'parameters_gnc.in';
filename_env  = 'parameters_env.json';

% Retrieve simulation time settings
[sim_start, sim_duration] = read_time_settings(filename_time); 

% Retrieve simulation environment settings
env = read_environment(filename_env, sim_start, sim_duration);

%% SPACESCRAFT Settings
% - Load from file the initial orbit parameters
filename_orbit = 'parameters_orbit_1_23042020.in';       

% Retrive orbital elements from file.
[~, coe] = read_orbit(filename_orbit);       

% Save the constant chaser COE
const_coe = reshape(coe(1:5), 5, 1);

% Compute Orbital Period (s) and Mean Motion (rad/s)
T_orbit = period(coe(1), env.GM_EARTH);
W_orb = mean_motion(coe(1), env.GM_EARTH);

% Set w threshold (rad/s)
w_threshold = 10*W_orb;

% Read TLEs
% tle_c = read_tle_lines('TLE_EXAMPLE.txt');

% Convert COE to a Cartesian state (position, velocity) in (m, m/s)
x_sc = coe2car(coe, env.GM_EARTH);

%% SIMULATION DATA

if n_sim == 1
    sim_type = 'OS';
    sim_mode = 'accelerator';

elseif n_sim > 1
    sim_type = 'MC';
    % rapid not working because of algebraic loops
    % assess possibility to remove algebraic loop with simplified algorithms 
    % in MC simulations
    % sim_mode='rapid';
    sim_mode = 'accelerator';

else
    error('Number of runs must be 1 or n>1')

end

% Initialize output vector
output = cell(n_sim, 1);

%% Set Single Simulation Parameters
tic

for i = 1:n_sim

    %% DKE Parameters Definition
    
    set_dke_parameters_sc(sim_type);
    set_sensors_parameters;
    set_actuators_parameters;
    
    %% GNC Parameters Definition 

    [TUN_NAV, NTUN_NAV, S_TIME_NAV] = NAV_parameters;
    [TUN_GUI, NTUN_GUI, S_TIME_GUI] = GUI_parameters;
    [TUN_CTR, NTUN_CTR, S_TIME_CTR] = CTR_parameters;
    [TUN_SENS_PP, NTUN_SENS_PP, S_TIME_SENS_PP] = SENS_PP_parameters;
    [TUN_ACTS_PP, NTUN_ACTS_PP, S_TIME_ACTS_PP] = ACTS_PP_parameters;
     
    %% BUS Parameters Definition 
    % Each function defines the busses that are associated with 
    % each specific subsystem.
    BUS_TEI_fun;
    BUS_DKE_fun;
    BUS_SENS_fun;
    BUS_ACTS_fun;
    BUS_OBC_fun;
    BUS_NAV_fun;
    BUS_GUI_fun;
    BUS_CTR_fun;
    BUS_SENS_PP_fun;
    BUS_ACTS_PP_fun;
    
    %% Simulation Initial conditions 
    switch sim_mode_type
        case 'detumbling'

            % Initial quaternion 
            sc.Q0 = rand(4, 1) - 0.5; % (-)
            sc.MaxW0_3sigma = 30; % (deg/s)
            % Initial angular velocity 
            sc.W0 = deg2rad(sc.MaxW0_3sigma/3)*randn(3, 1); % (rad/s)
            % Initial Wheel Angular Rates - 4 Wheels
            sc.WW0 = RWL_DKE.MaxWheelW*[0; 0; 0; 0]; % (rad/s)
            % True Anomaly
            sc.TA = 2*pi*rand(1); % (rad) - MC consider to put rand

        case 'pointing'
            % Initial quaternion
            sc.Q0 = rand(4, 1) - 0.5; % (-)
            % Initial angular velocity
            sc.MaxW0_3sigma = 10*W_orb; % With Margin 100% (rad/s) 
            sc.W0 = sc.MaxW0_3sigma/3*randn(3, 1);
            % Initial Wheel Angular Rates - 4 Wheels -- 50% of maximum
            sc.WW0 = 1*RWL_DKE.MaxWheelW*(rand(4, 1) - 0.5); % (rad/s)
            % True Anomaly
            sc.TA = coe(6); % (rad) - MC consider to put rand

        otherwise
            error('ADCS mode not valid')
    end
    
    %% DKE Initial Conditions 
    % -- Spacecraft Initial Conditions
    % Convert COE to an Inertial Cartesian State
    x_sc = coe2car([const_coe; sc.TA], env.GM_EARTH);
    sc.X0 = x_sc;
    sc.R0 = x_sc(1:3); % Position (m)
    sc.V0 = x_sc(4:6); % Velocity (m/s)

    % Random Initial Quaternion - Scalar Last
    sc.Q0 = sc.Q0./norm(sc.Q0);
    
 
    %% SIMULATION 
    % Add ALL sample times in simulation
    time_steps = [S_TIME_NAV.time_1, S_TIME_GUI.time_1, S_TIME_GUI.time_2, ...
                  S_TIME_CTR.time_1, S_TIME_CTR.time_2, S_TIME_ACTS_PP.time_1, ...
                  S_TIME_SENS_PP.time_1, MTQ_actuator_time, RWL_actuator_time , ...
                  THR_actuator_time];            
              
    sc.Q0 = rand(4, 1) - 0.5; % (-)
     
    % Initial angular velocity
    sc.MaxW0_3sigma = 10*W_orb; % With Margin 100% (rad/s) 
    sc.W0 = sc.MaxW0_3sigma/3*randn(3, 1);
    
    % Initial Wheel Angular Rates - 4 Wheels -- 50% of maximum
    sc.WW0 = 1*RWL_DKE.MaxWheelW*(rand(4, 1) - 0.5); % (rad/s)
    
    % True Anomaly
    sc.TA = coe(6); % (rad) - MC consider to put rand

    % Set Simulation Integration Step Time
    int_step_val = min(time_steps); 
    int_step = num2str(int_step_val);
    
    % Run Simulation
    if strcmp(sim_type, 'OS')
        fprintf(['Simulation initialized with:\n - Q0 = [%.2f, %.2f, %.2f, %.2f]\n', ...
            ' - W0 = [%.6f, %.6f, %.6f] deg/s\n - Ww0 = [%.2f, %.2f, %.2f, %.2f] rpm\n', ...
            ' - tsim = %.4f h = %i s\n - Repeatable = %i\n - Flex = %i\n' ...
            ' - ADS performance model disturbances = %i\n', ...
            ' - Actuator model disturbances = %i\n', ...
            ' - Solar panels: open = %i - BOL = %i\n', ...
            ' - sim_type = %s\n - gnc_mode = %i\n\n'], sc.Q0, (180/pi)*sc.W0,...
            (60/2/pi)*sc.WW0, sim_duration/3600, sim_duration, repeatable, flex_flag, ...
            gnc_pm_disturbance_flag, acts_disturbance_flag, sc.panels_open_flag, ...
            sc.bol_flag, sim_mode_type, gnc_mode_flag);
    else
        if i == 1
            fprintf(['MC Simulation initialized with:\n - n_sim = %i\n - tsim = %.4f', ...
                'h = %i s\n - Repeatable = %i\n - Flex = %i\n', ...
                ' - ADS performance model disturbances = %i\n', ...
                ' - actuator model disturbances = %i\n', ...
                ' - Solar panels: open = %i - BOL = %i\n', ...
                ' - sim_type = %s\n - gnc_mode = %i\n\n'], n_sim, ...
                sim_duration/3600, sim_duration, ...
                repeatable, flex_flag, gnc_pm_disturbance_flag, ...
                acts_disturbance_flag, sc.panels_open_flag, ...
                sc.bol_flag, sim_mode_type, gnc_mode_flag);
        end

        fprintf('MC Run number: %d out ot %d \n', i, n_sim);

    end
    
    simout = sim(simulator, 'simulationmode', sim_mode, 'FixedStep', int_step, ...
        'solver', 'ode5', 'StopTime', num2str(sim_duration, 16), ...
        'SrcWorkspace', 'current');
    
    %% Parse Output 
    % Data logging
    output{i}.r = get(simout, 'r');
    output{i}.q = get(simout, 'q');
    output{i}.w = get(simout, 'w');
    
    % ... add others as needed

end

clear simout

e_t_sim = toc;
fprintf('----- %d Simulations in %f s----\n ', n_sim, e_t_sim);

save(out_name, '-v7.3')
if n_sim > 1
    % Save also compressed output with 1/10 of points
    output_compressor(out_name, 0.1, 1)
end


