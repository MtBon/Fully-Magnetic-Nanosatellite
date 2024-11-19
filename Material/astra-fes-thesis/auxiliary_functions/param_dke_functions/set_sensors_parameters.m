function set_sensors_parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set Parameters for Sensors of ASTRA FES
%
% (c) ASTRA - Aerospace Science and Technology Dept. - PoliMi -
% 2020
%
% Needed parameters can be recovered from the calling function as
% PARAM=evalin('caller','PARAM');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------------------------------------------------------------%%
%%--------------------------- Sensors DATA --------------------------------------%%
%%-------------------------------------------------------------------------------%%
%% Sensor Sample Times
GYRO_sensor_time = 0.1; %10Hz

% ... add other sensors


%% Sensor Low Pass Frequencies - To model filters inside sensors and NOT in Post Processing (PP Filters to be set in GNC parameters)
% Cut off frequency, F, at -3dB is expressed in Hertz
% Cut off frequencies F<(Sensor_Frequency/2)
GYRO_lp_f = 2.5;
% ... add other sensors

%% Sensor Low Pass Butterworth Order
% Order of Butterworth Filter for Sensors
lp_order = 4;

%% Initialize Sensors Data
GYRO = HP_IMU_MS3025_init(GYRO_sensor_time); %#ok<NASGU>
% ... add other sensors

%% Sensor LP IIR Filters - These Filters are Embedded in Sensor Component
% lp_order_th order butterworth (lp_order as first input parameter)
% Cut off frequency, N, at -3dB is expressed in Hertz (N*xx_sensor_time*2 as second input parameter)
% Cut off frequency F<(Sensor_Frequency/2)
% Low pass filter ('low' as third input paramters)
[lp_zeros_GYRO,lp_poles_GYRO,lp_gain_GYRO] = butter(lp_order,GYRO_lp_f*GYRO_sensor_time*2,'low'); %#ok<ASGLU>
% ... add other sensors with defined internal LP filters

%% assign all variables to base workspace
myVarList=who;
for indVar = 1:length(myVarList)
    assignin('caller',myVarList{indVar},eval(myVarList{indVar}))
end

return
