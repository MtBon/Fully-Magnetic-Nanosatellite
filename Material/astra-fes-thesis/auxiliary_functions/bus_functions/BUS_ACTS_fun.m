function cellInfo = BUS_ACTS_fun(varargin)
% BUS_ACTS returns a cell array containing bus object information
%
% Optional Input: 'false' will suppress a call to Simulink.Bus.cellToObject
%                 when the MATLAB file is executed.
% The order of bus element attributes is as follows:
%   ElementName, Dimensions, DataType, SampleTime, Complexity, SamplingMode, DimensionsMode, Min, Max, DocUnits, Description

suppressObject = false;
if nargin == 1 && islogical(varargin{1}) && varargin{1} == false
    suppressObject = true;
elseif nargin > 1
    error('Invalid input argument(s) encountered');
end

%% Variables
% Call Variables Needed to Define Buses
MTQ_actuator_time=evalin('caller','MTQ_actuator_time');
RWL_actuator_time=evalin('caller','RWL_actuator_time');
THR_actuator_time=evalin('caller','THR_actuator_time');
%% Define Bus
cellInfo = { ...
    %% INPUT
    
    % No BUS in INPUT - FES EXAMPLE
    
    %% OUTPUT
    % ACTS BUS
    { ...
    'ACTS_OUT', ...
    '', ...
    '', ...
    'Auto', ...
    '-1', {...
    {'u', 3, 'double', THR_actuator_time, 'real', 'Sample', 'Fixed', [], [], '', ''}; ...
    {'d_mtq', 3, 'double', MTQ_actuator_time, 'real', 'Sample', 'Fixed', [], [], '', ''}; ...
    {'m_rwl', 4, 'double', RWL_actuator_time, 'real', 'Sample', 'Fixed', [], [], '', ''}; ...
    {'m_pert_ext', 3, 'double', RWL_actuator_time, 'real', 'Sample', 'Fixed', [], [], '', ''}; ...
    } ...
    } ...
    % ADD ADDITIONAL BUS RELATED TO THE BLOCK BELOW!!!
    }';
if ~suppressObject
    % Create bus objects in the MATLAB base workspace
    Simulink.Bus.cellToObject(cellInfo)
end
