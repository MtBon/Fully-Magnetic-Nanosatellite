function set_actuators_parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set Parameters for Actuators of ASTRA FES
%
% (c) ASTRA - Aerospace Science and Technology Dept. - PoliMi -
% 2020
%
% Needed parameters can be recovered from the calling function as
% PARAM=evalin('caller','PARAM');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
acts_disturbance_flag=evalin('caller','acts_disturbance_flag'); %#ok<*NASGU>
sim_type=evalin('caller','sim_type');

%% ------------------------------------------------------------------------------%%
%%--------------------------- Actuators DATA --------------------------------------%%
%%-------------------------------------------------------------------------------%%
%% Actuator Sample Times
MTQ_actuator_time=0.1; %10Hz - Typically the same of ACS sample time
RWL_actuator_time=0.1; %10Hz - Typically the same of ACS sample time
THR_actuator_time=0.1;   %10Hz  - TODO: tunami
% ... add other actuators


%% Initialize Actuators Data  
% ACT struct for data used in GNC - To be used in eventual Param GNC Function
% ACT_DKE struct for data used in DKE
[~,MTQ_DKE]=HP_MT_GST600_init(MTQ_actuator_time); %#ok<ASGLU> - are needed in the model of MTQ - if MTQ are present, MTQ_DKE shall be present 
[~,RWL_DKE]=HP_RW_GSW600_init(RWL_actuator_time); %#ok<ASGLU> - are needed in the dynamics and in the model of RWL - RWL_DKE shall be always present
% ... add other actuators if present

%% assign all variables to base workspace
myVarList=who;
for indVar = 1:length(myVarList)
    assignin('caller',myVarList{indVar},eval(myVarList{indVar}))
end

return
