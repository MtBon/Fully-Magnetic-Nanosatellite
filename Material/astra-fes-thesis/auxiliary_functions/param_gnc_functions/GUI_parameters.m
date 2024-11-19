function [TUN_STRUCT,NTUN_STRUCT,SAMPLE_TIME] = GUI_parameters
% It returns a struct array containing tunable struct parameters
% for GUI Block Models

%% MPC Parameters

% Tunable
TUN_STRUCT.Ts           = double(60);               % sample time
TUN_STRUCT.Ns           = double(100);              % finite horizon
TUN_STRUCT.S            = double(eye(6,6)*3e2);     % state weight matrix
TUN_STRUCT.R            = double(eye(3,3)*1e2);     % control weight matrix
TUN_STRUCT.uthreshold   = double(1e-3*ones(3,1));   % control threshold   


% Static - Non Tunable
TUN_STRUCT.Umax         = repmat(TUN_STRUCT.uthreshold,TUN_STRUCT.Ns,1);
TUN_STRUCT.Umin         = repmat(-TUN_STRUCT.uthreshold,TUN_STRUCT.Ns,1);
NTUN_STRUCT.W           = [TUN_STRUCT.Umax; -TUN_STRUCT.Umin]; % control constraint
NTUN_STRUCT.N           = [eye(3*TUN_STRUCT.Ns,3*TUN_STRUCT.Ns);-eye(3*TUN_STRUCT.Ns,3*TUN_STRUCT.Ns)]; % state constraint


%% Simulation parameters
NTUN_STRUCT.n       = evalin('caller', 'W_orb');
NTUN_STRUCT.re      = evalin('caller', 'env.RADIUS_EARTH');
NTUN_STRUCT.mu      = evalin('caller', 'env.GM_EARTH');
NTUN_STRUCT.j2      = evalin('caller', 'env.J2');

%% Reference parameters

% Tunable
% 1st hold orbit
TUN_STRUCT.ROEtar.Ho0   = single([0;0.0028127;0;0.0001406;0;0.0001406]);
TUN_STRUCT.holdTime.t0  = single(1e5);
% 1st inspection orbit
TUN_STRUCT.ROEtar.Io0   = single([2.576e-05;0.0028127;0;0.0001406;0;0.0001406]);

% Non-tunable
a                       = evalin('caller', 'coe(1)');
e                       = evalin('caller', 'coe(2)');
i                       = evalin('caller', 'coe(3)');
raan                    = evalin('caller', 'coe(4)');
aop                     = evalin('caller', 'coe(5)');
th                      = evalin('caller', 'coe(6)');
refOE                   = [1e3*a,e,i,raan,aop,th];




%% Sample time
SAMPLE_TIME.time_1=double(1); %[s]
SAMPLE_TIME.time_2=SAMPLE_TIME.time_1*TUN_STRUCT.Ts; %[s]
% ... add other sample times if present

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
