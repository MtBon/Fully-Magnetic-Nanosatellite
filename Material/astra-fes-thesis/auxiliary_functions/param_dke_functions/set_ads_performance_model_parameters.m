function set_ads_performance_model_parameters(sim_type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set Parameters for GNC Performance Model of ASTRA FES
%
% (c) ASTRA - Aerospace Science and Technology Dept. - PoliMi -
% 2018
%
% Requires the simulation type: OS or MC to set dispersion
% Others needed parameters can be recovered from the calling function as
% PARAM=evalin('caller','PARAM');
%
% V 2.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Flags
gnc_pm_disturbance_flag=evalin('caller','gnc_pm_disturbance_flag');

%% Parameters Set Up
if nargin == 0
    sim_type='OS'; %#ok<*NASGU>
end

% Set dispersions
switch sim_type
    case 'OS'
        errads=0.25/3; %25% - 3sigma
    case 'MC'
        errads=0.5/3; %50% - 3sigma
    otherwise
        error('ERROR: sim type must be ''OS'' or ''MC'' ')
end

%% ------------------------------------------------------------------------------%%
%%-------------------------- GNC Sample Time   ----------------------------------%%
%%-------------------------------------------------------------------------------%%
GNC_PM_Freq=5; %Hz
GNC_PM_sampletime=1/GNC_PM_Freq; %s
% ... add other sample times for performance models as needed


%% ------------------------------------------------------------------------------%%
%%----------------------  ADS Models Disturbances -------------------------------%%
%%-------------------------------------------------------------------------------%
%
% Can be used as EXAMPLES
%
if gnc_pm_disturbance_flag
    %ADS disturbances are ON   
    %Position Vector
    biasE=0.1*(sign(randn(3,1)).*ones(3,1)+errads*randn(3,1)); % 10cm - Within the inbox of the spacecraft - Bias in the GPS Performance Model
    stdE=1.5*(ones(3,1)+errads*randn(3,1));  % 1.5 m (1-sigma) - STD in the the GPS Performance Model
    seedE=floor(4925845*rand(3,1)); %Seed in the GPS Performance Model 
    
    %Velocity Vector
    biasV=10*(sign(randn(3,1)).*ones(3,1)+errads*randn(3,1)); % 10 m/s - Bias in the GPS Performance Model
    stdV=0.5*(ones(3,1)+errads*randn(3,1));  % 0.5 m/s (1-sigma) - STD in the the GPS Performance Model
    seedV=floor(56543*rand(3,1)); %Seed in the GPS Performance Model 
    
    %Sun Sensors
    biasSS=deg2rad(1)*(sign(randn(3,1)).*ones(3,1)+errads*randn(3,1)); %1 deg - Bias in the Sun Sensors Performance Model
    stdSS=deg2rad(2)*(ones(3,1)+errads*randn(3,1)); %2 deg - STD in the Sun Sensors Performance Model
    seedSS=floor(6778362*rand(3,1)); %Seed in the Sun Sensors Performance Model 
    
    %Magnetometers Sensors
    biasMag=1e-6*(sign(randn(3,1)).*ones(3,1)+errads*randn(3,1)); %1uT Bias in the Magnetometers Performance Model
    stdMag=5e-7*(ones(3,1)+errads*randn(3,1)); %500nT - STD in the Magnetometers Performance Model
    seedMag=floor(24518*rand(3,1)); %Seed in the Magnetometers Performance Model
    
    %Gyroscope
    BiasGyro=deg2rad(0.01)*(sign(randn(3,1)).*ones(3,1)+errads*randn(3,1)); %0.01 deg/s - Bias in the Gyro Performance Model (Bias Estimated)
    RRW=deg2rad(1e-6/3600)*(ones(3,1)+errads*randn(3,1)); %1e-6deg/h Variance RRW Noise Performance Model
    RRWnoiseSeed=floor(78874*rand(3,1)); %Seed RRW Noise Performance Model
    RRWnoiseTs=GNC_PM_sampletime; %RRW Noise Sample Time
    ARW=deg2rad(1e-6)*(ones(3,1)+errads*randn(3,1)); %1e-6 deg/S Variance ARW Noise Performance Model
    ARWnoiseSeed=floor(7610*rand(3,1)); %Seed ARW Noise Performance Model
    ARWnoiseTs=GNC_PM_sampletime; %ARW Noise Sample Time
    
    %ADS Performance Model
    biasAds=deg2rad(2)*(sign(randn(3,1)).*ones(3,1)+errads*randn(3,1)); %2 deg - Bias in the ADS Performance Model
    stdAds=deg2rad(0.5)*(ones(3,1)+errads*randn(3,1)); %0.5 deg - STD in the ADS Performance Model (Noise)
    seedAds=floor(245311*rand(3,1)); %Seed in the ADS Performance Model (Noise)
    stdAds_fluct=deg2rad(15/3600)*(ones(3,1)+errads*randn(3,1)); %15 deg/h - STD in the ADS Performance Model (Fluctuation)
    seedAds_fluct=floor(776637*rand(3,1)); %Seed in the ADS Performance Model (Noise)
    
    %Wheel Transducer Performance Model
    biasWw=5*2*pi/60*(sign(randn(4,1)).*ones(4,1)+errads*randn(4,1)); %5 rpm - Bias in the Wheel Transducer Performance Model
    stdWw=10*2*pi/60*(ones(4,1)+errads*randn(4,1)); %10 rpm - STD in the Wheel Transducer Performance Model
    seedWw=floor(10012*rand(4,1)); %Seed in the Wheel Transducer Performance Model 
    
    %Magnetorquer Temperature Sensor Performance Model
    biasTmtq=1.5*(sign(randn(1,1)).*ones(1,1)+errads*randn(1,1)); %1.5 degC - Bias in the Magnetorquer Temperature Sensor Performance Model
    stdTmtq=0.2*(ones(1,1)+errads*randn(1,1)); %0.2 degC - STD in the Magnetorquer Temperature Sensor Performance Model
    seedTmtq=floor(4890*rand(1,1)); %Seed in the Magnetorquer Temperature Sensor Performance Model
    
    %Magnetorquer PWM Current Sensor Performance Model
    biasImtq=1e-3*(sign(randn(1,1)).*ones(1,1)+errads*randn(1,1)); %1 mA (0.1% of Max PWM) - Bias in the Magnetorquer PWM Current Sensor Performance Model
    stdImtq=1e-2*(ones(1,1)+errads*randn(1,1)); %0.01 A (1% of Max PWM) - STD in the Magnetorquer PWM Current Sensor Performance Model
    seedImtq=floor(51190*rand(1,1)); %Seed in the Magnetorquer PWM Current Sensor Performance Model
else
    %ADS disturbances are OFF   
    %Position Vector
    biasE=0*ones(3,1); %Bias in the GPS Performance Model
    stdE=0*ones(3,1); %STD in the GPS Performance Model
    seedE=1*ones(3,1); %Seed in the GPS Performance Model
    
    %Velocity Vector
    biasV=0*ones(3,1); %Bias in the GPS Performance Model
    stdV=0*ones(3,1); %STD in the GPS Performance Model
    seedV=1*ones(3,1); %Seed in the GPS Performance Model
    
    %Sun Sensors
    biasSS=0*ones(3,1); %Bias in the Sun Sensors Performance Model
    stdSS=0*ones(3,1); %STD in the Sun Sensors Performance Model
    seedSS=1*ones(3,1); %Seed in the Sun Sensors Performance Model
    
    %Magnetometers Sensors
    biasMag=0*ones(3,1); %Bias in the Magnetometers Performance Model
    stdMag=0*ones(3,1); %STD in the Magnetometers Performance Model
    seedMag=1*ones(3,1); %Seed in the Magnetometers Performance Model
    
    %Gyroscope
    BiasGyro=0*ones(3,1); %Bias in the Gyro Performance Model
    RRW=0*ones(3,1); %Variance RRW Noise Performance Model
    RRWnoiseSeed=1*ones(3,1); %Seed RRW Noise Performance Model
    RRWnoiseTs=GNC_PM_sampletime; %RRW Noise Sample Performance Model
    ARW=0*ones(3,1); %Variance ARW Noise Performance Model
    ARWnoiseSeed=1*ones(3,1); %Seed ARW Noise Performance Model
    ARWnoiseTs=GNC_PM_sampletime; %ARW Noise Sample Time
    
    %ADS Performance Model
    biasAds=0*ones(3,1); %Bias in the ADS Performance Model
    stdAds=0*ones(3,1); %STD in the ADS Performance Model
    seedAds=1*ones(3,1); %Seed in the ADS Performance Model
    stdAds_fluct=0*ones(3,1); %STD in the ADS Performance Model (Fluctuation)
    seedAds_fluct=1*ones(3,1); %Seed in the ADS Performance Model (Noise)
    
    %Wheel Transducer Performance Model
    biasWw=0*ones(4,1); %Bias in the Wheel Transducer Performance Model
    stdWw=0*ones(4,1); %STD in the Wheel Transducer Performance Model
    seedWw=1*ones(4,1); %Seed in the Wheel Transducer Performance Model
    
    %Magnetorquer Temperature Sensor Performance Model
    biasTmtq=0*ones(1,1); %Bias in the Magnetorquer Temperature Sensor Performance Model
    stdTmtq=0*ones(1,1); %STD in the Magnetorquer Temperature Sensor Performance Model
    seedTmtq=1*ones(1,1); %Seed in the Magnetorquer Temperature Sensor Performance Model
    
    %Magnetorquer PWM Current Sensor Performance Model
    biasImtq=0*ones(1,1); %Bias in the Magnetorquer PWM Current Sensor Performance Model
    stdImtq=0*ones(1,1); %STD in the Magnetorquer PWM Current Sensor Performance Model
    seedImtq=1*ones(1,1); %Seed in the Magnetorquer PWM Current Sensor Performance Model
end

%% assign all variables to base workspace
myVarList=who;
for indVar = 1:length(myVarList)
    assignin('caller',myVarList{indVar},eval(myVarList{indVar}))
end

return
