function cellInfo = BUS_GUI_fun(varargin)

    % BUS_GUI returns a cell array containing bus object information
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
    S_TIME_GUI = evalin('caller', 'S_TIME_GUI');

    %% Define Bus
    cellInfo = { ...
        %% INPUT
        % GUI_NAV_IN BUS
        { ...
        'GUI_NAV_IN', ...
        '', ...
        '', ...
        'Auto', ...
        '-1', {...
        {'x_est', 6, 'single', S_TIME_GUI.time_1, 'real', 'Sample', 'Fixed', [], [], '', ''}; ...
        } ...
        } ...
        
        % MODE MANAGER
        { ...
        'GUI_MM_IN', ...
        '', ...
        '', ...
        'Auto', ...
        '-1', {...
        {'dig_clock', 1, 'double', S_TIME_GUI.time_1, 'real', 'Sample', 'Fixed', [], [], '', ''}; ...
        } ...
        } ...

        %% OUTPUT
        % GUI MSG
        { ...
        'GUI_MSG_OUT', ...
        '', ...
        '', ...
        'Auto', ...
        '-1', {...
        {'gui_flag', 1, 'uint8', S_TIME_GUI.time_1, 'real', 'Sample', 'Fixed', [], [], '', ''}; ...
        } ...
        } ...
        
        % GUI DATA
        { ...
        'GUI_DATA_OUT', ...
        '', ...
        '', ...
        'Auto', ...
        '-1', {...
        {'uMPC', [3,1], 'single', S_TIME_GUI.time_2, 'real', 'Sample', 'Fixed', [], [], '', ''}; ...
        } ...
        } ...
        
        % ADD ADDITIONAL BUS RELATED TO THE BLOCK BELOW!!!
        }';
    if ~suppressObject
        % Create bus objects in the MATLAB base workspace
        Simulink.Bus.cellToObject(cellInfo)
    end

end