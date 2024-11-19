function [RWL,RWL_DKE]=HP_RW_GSW600_init(actuator_time)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RW errors initialization %
% GSW-600
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
acts_disturbance_flag=evalin('caller','acts_disturbance_flag');
sim_type=evalin('caller','sim_type');

% Set dispersions
switch sim_type
    case 'OS'
        erract=0.25/3; %25% - 3sigma
    case 'MC'
        erract=0.5/3; %50% - 3sigma
    otherwise
        error('ERROR: sim type must be ''OS'' or ''MC'' ')
end

%% Physical Parameters
% Single Wheel Inertia - kg*m2
% Assuming a radius of the wheel equal to r=20mm (cfr. DS at p.11)
% the mass of the wheel is ~2*I/r^2=0.151 kg
I=3.0239e-05; %19mNms at 6000rpm 

% Wheels Mounting Matrix - Pyramid at 60deg with respect horizontal
A=0.5*[-sqrt(3), sqrt(3), 0, 0; 0, 0, -sqrt(3), sqrt(3); 1, 1, 1, 1,];

% Wheel Friction
if acts_disturbance_flag
%Wheels Dead Zone
DeadZone=150*(2*pi/60); %[rad/s] - Stribeck Model
ZeroCross=10*(2*pi/60); %[rad/s] - Karnopp Model
Stiction=1.75; %[%/100] - Stiction Model
% Friction
TorqueBias=4e-5; % Coulomb Model term
TorqueFriction=3e-6; %Viscous Model term
else
%Wheels Dead Zone
DeadZone=1; %[rad/s] - Stribeck Model
ZeroCross=1; %[rad/s] - Karnopp Model
Stiction=0; %[%/100] - Stiction Model
% Friction
TorqueBias=0; % Coulomb Model term
TorqueFriction=0; %Viscous Model term   
end

% Wheel Limits
MaxTorque=1.5e-3; %[Nm] - Maximum Wheel Control Torque
MaxWheelW=2*pi*6000/60; %[rad/s] - Maximum Wheel Speed

%% Parameters for ACS
%%%%% RWL INERTIA PARAMETERS
% Wheels Inertia - kg*m2
RWL.I=single(I);
% Wheels Inertia Matrix - kg*m2
RWL.AI=single(A*I);

%%%%% RWL PERFORMANCE PARAMETERS
% Wheel  Torque Box Limits - Nm - rad/s
RWL.MaxTorque=single(MaxTorque);
RWL.MaxTorque_MaxW=single(MaxTorque/3);
RWL.MaxWheelW=single(MaxWheelW); 
RWL.MaxWheelW_MaxT=single(2*pi*3500/60); 

%%%%% RWL ACTUATION
% Wheels Cycle Time - S
RWL.Ts=single(5*actuator_time); %[~5/10 times the ACS loop time]

%%%%% RWL DISTRIBUTION MATRICES
% Wheels Distribution - Pseudo Inverse
RWL.invA=pinv(A); % All RWL Working
RWL.invA(abs(RWL.invA)<1e-12)=0;
RWL.invA=single(RWL.invA);
% Fail RWL_1 (-x)
A_sF=A;
A_sF(:,1)=zeros(3,1);
RWL.invA_F1=pinv(A_sF);
RWL.invA_F1(abs(RWL.invA_F1)<1e-12)=0;
RWL.invA_F1=single(RWL.invA_F1);
% Fail RWL_2 (+x)
A_sF=A;
A_sF(:,2)=zeros(3,1);
RWL.invA_F2=pinv(A_sF);
RWL.invA_F2(abs(RWL.invA_F2)<1e-12)=0;
RWL.invA_F2=single(RWL.invA_F2);
% Fail RWL_3 (-y)
A_sF=A;
A_sF(:,3)=zeros(3,1);
RWL.invA_F3=pinv(A_sF);
RWL.invA_F3(abs(RWL.invA_F3)<1e-12)=0;
RWL.invA_F3=single(RWL.invA_F3);
% Fail RWL_4 (+y)
A_sF=A;
A_sF(:,4)=zeros(3,1);
RWL.invA_F4=pinv(A_sF);
RWL.invA_F4(abs(RWL.invA_F4)<1e-12)=0;
RWL.invA_F4=single(RWL.invA_F4);

%%%%% RWL SATURATION DETCTION
% Single RWL Saturation Detection
RWL.S_DESAT_LIMIT=single(MaxWheelW*0.95); %[rad/s]
% Total RWL Saturation Detection - Limits and Delay Time
RWL.T_DESAT_LIMIT=single(MaxWheelW*0.75); %[rad/s]
RWL.T_DESAT_DELAY=single(10); %[s]

% %%%%% RWL FRICTION COMPENSATION MODEL
% RWL.DeadZone=DeadZone;
% RWL.ZeroCross=ZeroCross;
% RWL.Stiction=Stiction;
% RWL.TorqueBias=TorqueBias;
% RWL.TorqueFriction=TorqueFriction;

%% Parameters for DKE
% Disturbances
if acts_disturbance_flag
    %Hardware disturbances are ON
    % Wheels Dead Zone
    RWL_DKE.DeadZone=DeadZone*(ones(4,1)+erract*randn(4,1)); %[rad/s] - W_s: Stribeck Wheel Speed ~W_DZ/2 (N.B. RWL Dead Zone (300 RPM))
    RWL_DKE.ZeroCross=ZeroCross*(ones(4,1)+erract*randn(4,1)); %[rad/s] - W_zc: Velocity where Static Friction has Sign Change (Karnopp Model)
    RWL_DKE.Stiction=Stiction*(ones(4,1)+erract*randn(4,1)); %[%/100] - Stiction Gain Compared to Coulomb Force (i.e. Fs/Fc-1)
    % Friction
    RWL_DKE.TorqueFrictionGain=2.5*(ones(4,1)+erract*randn(4,1)); % Gain to increase friction in case of RWL failure
    RWL_DKE.TorqueBias=TorqueBias*(ones(4,1)+erract*randn(4,1)); %Bias in the Wheels Torque - Coulomb term
    RWL_DKE.TorqueFriction=TorqueFriction*(ones(4,1)+erract*randn(4,1)); %Bias in the Wheels Torque - Viscous term    
    % Noise/Imbalance - ISO1940-1 to Quality G 0.4 at 5000 rpm
    % G0.4 --> e*W=0.4 mm/s --> Static unbalance: U=1000*(e*W)*m/W 
    % W is used at quality level definition: W=5000rpm=523.59rad/s
    % U=1000*0.4*0.151/523.6=0.115 g*mm ~ 1.15e-7 kg*m
    % Assuming from literature a Ud~75*Us --> Ud=8.625g*mm2 ~ 8.625e-9 kg*m2 
    % Noise=Ud*Ww^2 --> MAX N_XY=3.4e-3 Nm (X and Y)
    % Assuming from literature MAX N_Z=0.025*N_XY=8.5e-5 Nm (Z)
    % Noise/Imbalance
    RWL_DKE.NoiseStd=1/3*(0.025*8.625e-9)*(ones(4,1)+erract*randn(4,1)); %Ud - Imbalance Noise Coefficient
    RWL_DKE.NoiseSeed=floor(666*rand(4,1)); %Seed in the Wheels Torque
    % Mounting Error
    ax_err=randn(1,3);
    RWL_DKE.MountErr = vrrotvec2mat([ax_err./norm(ax_err),deg2rad(0.1*randn)]); %0.1 deg - Mounting error matrix
else
    %Hardware disturbances are OFF
    %Wheels Dead Zone
    RWL_DKE.DeadZone=DeadZone*ones(4,1); %W_s: Stribeck Wheel Speed
    RWL_DKE.ZeroCross=ZeroCross*ones(4,1); %W_zc: Velocity where Static Friction has Sign Change
    RWL_DKE.Stiction=Stiction*ones(4,1); %Stiction Gain Compared to Coulomb Force
    % Friction
    RWL_DKE.TorqueFrictionGain=0*ones(4,1); % Gain to increase friction in case of RWL failure
    RWL_DKE.TorqueBias=TorqueBias*ones(4,1); %Bias in the Wheels Torque
    RWL_DKE.TorqueFriction=TorqueFriction*ones(4,1); %Bias in the Wheels Torque
    % Noise/Imbalance
    RWL_DKE.NoiseStd=0*ones(4,1);  %Ud - Imbalance Noise Coefficient
    RWL_DKE.NoiseSeed=1*ones(4,1); %Seed in the Wheels Torque
    % Mounting Error
    RWL_DKE.MountErr=eye(3); %Mounting error matrix
end

%Wheels Speed Control Accuracy
RWL_DKE.AccuracyW=2*pi*0.5/60; %[rad/s] - Wheel Speed Control Accuracy

%Wheels Speed Control Gain
RWL_DKE.KPT_MCU=1.5; %MCU Proportonal Gain Torque - Simplified MCU Control
RWL_DKE.KIT_MCU=7.5; %MCU Integral Gain Torque - Simplified MCU Control 
RWL_DKE.KPW_MCU=1e-6; %MCU Proportonal Gain W - Simplified MCU Control 
RWL_DKE.KIW_MCU=1e-9; %MCU Integral Gain W - Simplified MCU Control 
RWL_DKE.AWU_MCU=1; %MCU Anti Wind-up

%Wheels Speed/Torque Limits
RWL_DKE.MaxTorque=2.5e-3;%[Nm] RWL MCU MaxTorque
RWL_DKE.MaxTorque_MaxW=2.0e-3; %[Nm] -  RWL MCU MaxTorque at Practical Maximum speed
RWL_DKE.MaxTorque_Out=1.55e-3;%[Nm] RWL MaxTorque Output
RWL_DKE.MaxWheelW=MaxWheelW; %[rad/s] - Maximum Wheel Speed
RWL_DKE.MaxWheelW_Trq=2*pi*5900/60; %[rad/s] - Maximum Wheel Speed for Practical Torque Generation (6000-100 rpm)
RWL_DKE.MaxWheelW_Back=2*pi*5800/60; %[rad/s] - Maximum Wheel Speed for backward velocity (6000-200 rpm)
RWL_DKE.MaxWheelW_Tlin=2*pi*3500/60; %[rad/s] - Maximum Wheel Speed for Torque Linearity

%Wheels Power Supply Model
RWL_DKE.V=1.1*5; %[V] Wheel Power Supply: Voltage (5V and 10% margin)
RWL_DKE.Imin=0.03; %[A] Wheel Power Supply: Minimum Current (0.03 A=30mA)

% Same Quantities Described Above
RWL_DKE.I=I;
RWL_DKE.AI=A*I;

