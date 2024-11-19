

function data = read_json(filename)

    % Load a JSON file and extract its data to a MATLAB structure.
    %
    % Parameters
    % ----------
    %   filename: string
    %       Target file name.
    %   
    % Returns
    % -------
    %   data: struct
    %       JSON data as a MATLAB structure.

    filetext = fileread(filename);
    data = jsondecode(filetext);

end