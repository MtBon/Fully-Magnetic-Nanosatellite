function cellInfo = BUS_SENS_fun(varargin)

    % BUS_SENS returns a cell array containing bus object information
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
    GYRO_sensor_time = evalin('caller', 'GYRO_sensor_time');
    
    %% Define Bus
    cellInfo = { ...
        %% INPUT

        % No BUS in INPUT

        %% OUTPUT
        % SENS BUS
        { ...
        'SENS_OUT', ...
        '', ...
        '', ...
        'Auto', ...
        '-1', {...
        {'w_gyro', 3, 'single', GYRO_sensor_time, 'real', 'Sample', 'Fixed', [], [], '', ''}; ...
        {'t_gyro', 1, 'single', GYRO_sensor_time, 'real', 'Sample', 'Fixed', [], [], '', ''}; ...
        } ...
        } ...
        % ADD ADDITIONAL BUS RELATED TO THE BLOCK BELOW!!!
        }';
    
    if ~suppressObject
        % Create bus objects in the MATLAB base workspace
        Simulink.Bus.cellToObject(cellInfo)
    end
    
end
