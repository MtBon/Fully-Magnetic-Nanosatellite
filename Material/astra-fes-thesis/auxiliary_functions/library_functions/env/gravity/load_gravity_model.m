
%{
    load_gravity_model(filename, maxdeg)

Parse an official ICGEM gravity model file and extract the data
requested for the computation of the spherical harmonic expansion. 

# Inputs 
- `filename`: *string* -- Gravity model filename. 
- `maxdeg`: *integer* -- Maximum desired expansion degree. The coefficients
    associated to higher orders will not be extracted from the file. If the 
    specified degree is higher than that available in the file, an error is 
    thrown. 

# Outputs 
- `GRAV`: *struct* -- Data structure storing the following properties: 
    - *maxdeg*: **integer** -- Maximum expansion degree.
    - *GM*: **double** - Model gravity constant, in km³/s².
    - *R*: **double** - Model reference radius, in km.
    - *Chn*: **double (N, 1) ** - Normalised cosine coefficients.
    - *Shn*: **double (N, 1) ** - Normalised sine coefficients.
    - *K*: **double (M, 1) ** - Pre-computed scale factor coefficients.
    
# References
- [ICGEM](http://icgem.gfz-potsdam.de/ICGEM-Format-2011.pdf)

%}

function GRAV = load_gravity_model(filename, maxdeg)

    % Parse an official ICGEM gravity model file and extract the data
    % requested for the computation of the spherical harmonic expansion. 
    %
    % Parameters
    % ----------
    %   filename: string 
    %       Gravity model filename
    %   maxdeg: integer
    %       Maximum desired expansion degree. The coefficients associated 
    %       to higher orders will not be extracted from the file. If the 
    %       specified degree is higher than that available in the file, 
    %       an error is thrown. 
    %
    % Returns 
    % -------
    %   GRAV: struct 
    %       Data structure storing the following properties:
    %           - maxdeg (integer): Maximum expansion degree.
    %           - GM: double - Model gravity constant, in m³/s².
    %           - R: double - Model reference radius, in m.
    %           - Chn: double (N, 1) - Normalised cosine coefficients.
    %           - Shn: double (N, 1) - Normalised sine coefficients.
    %           - K: double (M, 1) - Pre-computed scale factor coefficients.
    % 
    % References 
    % ----------
    %   [1]: ICGEM, http://icgem.gfz-potsdam.de/ICGEM-Format-2011.pdf

    %#codegen 

    fid = fopen(filename, 'r');
    
    % Parse the file header and retrieve the gravity model properties
    [header, header_line_end] = parse_header(fid);
    [GRAV.GM, GRAV.R, GRAV.maxdeg, err] = parse_header_data(header);
    
    % Check that the model has enough coefficients for the desired degree
    if GRAV.maxdeg < maxdeg
        error(['The ICGEM model does not contain enough coefficients ' ...
               'for the desired %s degree'], maxdeg);
    end

    % Update with maximum desired degree and associated vector size
    GRAV.maxdeg = maxdeg;
    vsize = 0.5*(maxdeg+1)*(maxdeg+2);

    Chn = zeros(vsize, 1);
    Shn = zeros(vsize, 1);

    cline = header_line_end + 1;
    is_finished = false; 

    % Parse all the coefficients
    while ~feof(fid) && ~is_finished
        line = strsplit(fgetl(fid), ' ');
       
        % Check and update the number of columns 
        if cline == header_line_end + 1
            ncols = length(line);
            validate_columns(err, ncols);
        end

        % We only support the gfc coefficients.
        if ~strcmp(line(1), "gfc")
            continue;
        end

        % Retrieve current degree and order
        l = str2double(replace_exp(line(2)));
        m = str2double(replace_exp(line(3)));    
        idx = l*(l+1)/2 + m + 1;
        
        % Extract the coefficients
        Chn(idx) = str2double(replace_exp(line(4)));
        Shn(idx) = str2double(replace_exp(line(5)));

        is_finished = idx >= vsize;
        
        cline = cline + 1;
    end
    
    % Check to ensure that the number of columns is as expected.
    GRAV.Chn = Chn;
    GRAV.Shn = Shn;
    
    % Pre-compute scale factors which hold irrespective of the position
    GRAV.K = scale_factors(maxdeg);

end


function K = scale_factors(maxdeg)

    % Compute the scale factors up to the desired degree and
    % efficiently store them in vector form.

    %#codegen 

    vsize = 0.5*(maxdeg+3)*(maxdeg+4);

    K = zeros(vsize, 1);
    K(2) = 1; 

    for j = 2:maxdeg+2
        IDx = j + 1; 
        crt = 0.5*IDx*(IDx + 1);

        for k = 0:j-1
            if k == 0
                K(crt-j+k) = sqrt(0.5*(j+k+1)*(j-k));
            else 
                K(crt-j+k) = sqrt((j+k+1)*(j-k));
            end
        end
    end

end


function [header, cline] = parse_header(fid)

    % Extract from a ICGEM file the header string, which should be 
    % between the "begin_of_head" and "end_of_head" keywords.

    % Since the begin_of_head keyword is optional to maintain
    % compatibility with older releases, we initially push inside 
    % header all the lines, which are then removed if the keyword is 
    % found

    cline = 0; 
    
    head_end_found = false;

    % We do not pre-allocate this array, because this is not 
    % a time-critical function.
    header = "";
    
    % Parse header and look for the start and end keywords
    while ~feof(fid)
        cline = cline + 1;
        line = strsplit(fgetl(fid), ' ');

        % Empty lines 
        if length(line) < 1
            continue;
        end

        if strcmp(line(1), "begin_of_head")
            header = ""; 
            continue;
        end

        if strcmp(line(1), "end_of_head")
            head_end_found = true;
            break;
        end

        % Push current expression into the header
        header = strcat(header, join(line, " "), " ");

    end

    % End of head keyword is mandatory! 
    if ~head_end_found
        error("Invalid ICGEM format: keyword `end_of_head` was not found!");
    end

end


function [GM, R, maxdeg, err] = parse_header_data(header)

    % Extract from a header string the gravity model properties.

    %#codegen 

    % Keywords for data to be looked for
    keywords = ["product_type", "gravity_constant", "radius", ...
                "max_degree", "errors"];
    
    for k = 1:length(keywords)
        m = regexp(header, sprintf("%s .*[\\d\\w]", keywords(k)), 'match');

        if isempty(m)
            error("Invalid ICGEM format: keyword %s not found!", keywords(k));
        end

        matched = strsplit(m, " ");

        switch keywords(k) 
            case "product_type"
                isgravmodel = strcmp(matched(2), "gravity_field");

                if ~isgravmodel
                    error("Only ICGEM gravity_field products are supported");
                end

            case "gravity_constant"
                % Gravity constant in m³/s²
                GM = str2double(replace_exp(matched(2)));
                
            case "radius"
                % Radius in meters!
                R = str2double(replace_exp(matched(2)));

            case "max_degree"
                maxdeg = str2double(replace_exp(matched(2)));
                
            case "errors"
                err = matched(2);

                if ~any(err == ["no", "calibrated", "formal", ...
                                    "calibrated_and_formal"])

                    error("Invalid ICGEM error type!")
                end
        end

    end

end


function validate_columns(errtype, ncols)
    
    % Validate the number of columns found according to the specified error
    % type. Throw an error if a mismatch happens.

    %#codegen

    if strcmp(errtype, "no") 
        if ncols ~= 5
            error(['Invalid ICGME format: the coefficients table must have ' ...
                '5 columns when the keyword error is %s'], errtype);
        end

    elseif (strcmp(errtype, "calibrated") || strcmp(errtype, "formal")) 
        if ncols ~= 7
            error(['Invalid ICGME format: the coefficients table must have ' ...
               '7 columns when the keyword error is %s'], errtype);
        end
        
    elseif ncols ~= 9
        error(['Invalid ICGME format: the coefficients table must have ' ...
               '9 columns when the keyword error is %s'], errtype);
    end

end

function sub = replace_exp(raw)
    % Replace the D and d in the string `raw`
    % with the exponential `E`

    sub = strrep(strrep(raw, 'D', 'E'), 'd', 'E');

end
