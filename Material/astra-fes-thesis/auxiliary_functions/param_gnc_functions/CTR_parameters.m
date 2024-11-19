function [TUN_STRUCT,NTUN_STRUCT,SAMPLE_TIME] = CTR_parameters
% it returns a struct array containing tunable struct parameters
% for CTR Block Models

% Tunable
TUN_STRUCT.uthreshold   = single(1e-2*ones(3,1));
TUN_STRUCT.thrustTime   = single(60);
TUN_STRUCT.thrustRest   = double(600);
TUN_STRUCT.Ts           = double(60);               % MPC sample time
TUN_STRUCT.u0           = single(zeros(3, 1));
TUN_STRUCT.maxAngle     = single(20); % [deg]

% Sample time
SAMPLE_TIME.time_1=1; %[s]
SAMPLE_TIME.time_2=SAMPLE_TIME.time_1*TUN_STRUCT.Ts; %[s]
% ... add other sample times if present


% Static - Non Tunable - LQR
% Clohessy-Wiltshire linear state space model
n       = evalin('caller', 'W_orb');
A       = [zeros(3,3) eye(3,3);
            3*n^2 0 0 0 2*n 0;
            0 0 0 -2*n 0 0;
            0 0 -n^2 0 0 0];
B       = [zeros(3,3);eye(3,3)];
C       = eye(6,6);
D       = zeros(6,3);

% Continous state space system
cwsys   = ss(A,B,C,D);
cwsys_d = c2d(cwsys,SAMPLE_TIME.time_1,'tustin');

TUN_STRUCT.S        = double(eye(6,6)*3e2);
TUN_STRUCT.R        = double(eye(3,3)*1e2);

NTUN_STRUCT.Klqr    = lqr(cwsys_d,TUN_STRUCT.S,TUN_STRUCT.R);


% ....

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
