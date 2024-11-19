function cellInfo = BUS_OBC_fun(varargin)
% BUS_OBC returns a cell array containing bus object information
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

%% Define Bus
cellInfo = { ...
    %% INPUT
    
    % No BUS in INPUT - FES EXAMPLE
    
    %% OUTPUT
    % OBC 1 BUS
    { ...
    'OBC_1_OUT', ...
    '', ...
    '', ...
    'Auto', ...
    '-1', {...
    {'mode', 1, 'uint16', -1, 'real', 'Sample', 'Fixed', [], [], '', ''}; ...
    {'opt_1', 1, 'boolean', -1, 'real', 'Sample', 'Fixed', [], [], '', ''}; ...
    {'dig_clock', 1, 'double', -1, 'real', 'Sample', 'Fixed', [], [], '', ''}; ...
    } ...
    } ...
    % ADD ADDITIONAL BUS RELATED TO THE BLOCK BELOW!!!
    }';
if ~suppressObject
    % Create bus objects in the MATLAB base workspace
    Simulink.Bus.cellToObject(cellInfo)
end
