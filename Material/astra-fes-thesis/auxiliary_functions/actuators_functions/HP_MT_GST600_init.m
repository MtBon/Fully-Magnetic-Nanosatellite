function [MTQ,MTQ_DKE]=HP_MT_GST600_init(~)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MTORQ errors initialization % 
% GST-600
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

%% Parameters for PWM
PWM_V=3.3; %V - Supply Voltage

%% Actuators Magnetic Charachteristics
MAG_Range=[0.31;0.31;0.34]; %[Am2] - Magnetic Dipole at 3.3V and 20degC

%% Actuator Electrical Charachteristics
COIL_Resistance=[28.5;28.5;25.5]; %[Ohm] - Electrical Resistance at 20 degC
COIL_T_ref=20; %[degC] - Reference temperature (20 degC)
COIL_T_coeff=[0.112;0.112;0.104]; %[Ohm/degC] - Resistance Temperature Coefficient
% COIL_Max_Resistance=COIL_Resistance+COIL_T_coeff*(100-COIL_T_ref); %[Ohm] - Electrical Resistance at 100 degC
% COIL_Min_Resistance=COIL_Resistance+COIL_T_coeff*(-50-COIL_T_ref); %[Ohm] - Electrical Resistance at -50 degC

%% Parameters for ACS
% Actuator Performances Charachteristics
MTQ.range=single(MAG_Range); %Magnetic Dipole Range
MTQ.V_range=single(PWM_V); %[V] - Maximum PWM Voltage

% Actuator Resistance Estimation Parameters
MTQ.R_0=COIL_Resistance; %Initial Resistance Delay
MTQ.T_0=COIL_T_ref;
MTQ.T_coeff=COIL_T_coeff;
% MTQ.Mean_R_Samples=100; %Number of Samples to Average R
% MTQ.Max_R=single(COIL_Max_Resistance);
% MTQ.Min_R=single(COIL_Min_Resistance);

% Actuator Scaling Coefficiet - nA
MTQ.nA_coeff=single((MAG_Range.*COIL_Resistance)/PWM_V);

%% Parameters for DKE
% Temperature
MTQ_DKE.T=40*(ones(1,1)+erract*randn(1,1)); %[degC] Typical operating temperature - from TCS
MTQ_DKE.dT=20*(ones(1,1)+erract*randn(1,1)); %[degC] Typical temperature fluctuations - from TCS

% Power
MTQ_DKE.P_marg=1.1; %[%/100] 10% margin 

% PWM
MTQ_DKE.PWM_V=PWM_V;
MTQ_DKE.PWM_RES = PWM_V/2^16; % 16-Bit - Resolution
%
% ADDITIONAL PWM DATA
%
% MTQ_DKE.PWM_I=1; %A - Total and Channel Current Output
% MTQ_DKE.PWM_f=250e3; %[Hz] - Frequency
% Above quantities unused for practical applications:
% - PWM max current output of 1A is higher than maximum MTQ current draw 
%   (i.e. ~25 Ohm at 3.3 --> 0.13A per channel --> 0.4A per all 3 channels)
% - PWM frequency is extremely high and not influent since f_ACS is much lower

% Disturbances
if acts_disturbance_flag
    %Hardware disturbances are ON
    %Magnetorquers - PWM
    MTQ_DKE.PWM_noise=(PWM_V/1000)*(ones(3,1)+erract*randn(3,1)); %0.1% Max - STD in the PWM Noise
    MTQ_DKE.PWM_bias=(PWM_V/2000)*(sign(randn(3,1)).*ones(3,1)+erract*randn(3,1)); %0.05% Max - Bias in the PWM
    MTQ_DKE.PWM_noiseseed=floor(4978845*rand(3,1)); %Seed in the PWM Noise
    %Magnetorquer - COILS
    MTQ_DKE.residual=1e-3*(ones(3,1)+erract*randn(3,1)); %[Am2] 1mAm2 - Residual Dipole in the Magnetorquer  
    MTQ_DKE.sfn=(10000*1e-6)*(ones(3,1)+erract*randn(3,1))./(MTQ.range); %[ppm] 10000ppm - Scale factor nonlinearity in the Magnetorquer Dipole   
else
    %Hardware disturbances are OFF
    %Magnetorquers - PWM
    MTQ_DKE.PWM_noise=0.0*ones(3,1); %STD in the PWM Noise
    MTQ_DKE.PWM_bias=0.0*ones(3,1); %Bias in the PWM
    MTQ_DKE.PWM_noiseseed=1*ones(3,1); %Seed in the PWM Noise
    %Magnetorquer - COILS
    MTQ_DKE.residual=0.0*ones(3,1); %Residual Dipole in the Magnetorquer   
    MTQ_DKE.sfn=0.0*ones(3,1); %Scale factor nonlinearity in the Magnetorquer Dipole  
end

% Same Quantities Described Above
MTQ_DKE.resistance=COIL_Resistance;
MTQ_DKE.T_ref=COIL_T_ref;
MTQ_DKE.T_coeff=COIL_T_coeff;
MTQ_DKE.nA_coeff=MTQ.nA_coeff;