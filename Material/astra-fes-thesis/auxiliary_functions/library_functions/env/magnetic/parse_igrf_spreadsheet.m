function parse_igrf_spreadsheet(filename, row_start, col_start, output_name)

    % Parse the official IGRF excel spreadsheet containing the model coefficients and 
    % export those corresponding to the desired 5-year period to a text file. This 
    % function assumes that the coefficients in the excel file are Schmidt semi-normalised
    % and automatically converts them to the Gaussian normalisation. 
    % An output file will be created containing the cosine and sine coefficients. The columns 
    % of the file are organised as follows: [degree, order, g, dg, h, dh]
    %
    % Parameters
    % ----------
    %  filename: string 
    %         Name of the input spreadsheet file.
    %  row_start: integer
    %         The number of the first row that has valid coefficients. 
    %  col_start: integer
    %         The number of the column containing the reference epoch coefficients. If it 
    %         equals the second to last column, the function assumes the last column has 
    %         the secular variation coefficients. Otherwise, the secular variation is 
    %         obtained from a linear interpolation of the coefficients of the current and 
    %         next epochs.
    % 
    % Notes
    % -----
    % This function is only intended to be used to generate the new reference text file.
    % This file is then properly parsed by the mag_model function.
    % 
    % Example
    % -------
    %   parse_igrf_spreadsheet('IGRF13coeffs.xls', 5, 28, 'IGRF20');
    % 
    % References
    % ----------
    % [1] J. R. Wertz, Spacecraft Attitude Determination and Control, 1978, Appendix H.

    if row_start < 1 
        error("Starting row must be >= 1.");
    end

    if col_start < 1 
        error("Starting column must be >= 1.");
    end

    opts = detectImportOptions(filename);
    ncols = length(opts.VariableOptions); 
    
    if col_start >= ncols
        error("Starting column must be < than the number of available columns in the file.");
    end

    opts.DataRange = sprintf('A%d', row_start);

    %  Extract the Schmidt semi-normalised spherical harmonics coefficients
    mag_tab = readtable(filename, opts); 
    
    letter = table2array(mag_tab(:, 1));
    deg = table2array(mag_tab(:, 2)); 
    ord = table2array(mag_tab(:, 3)); 

    coeff = table2array(mag_tab(:, col_start));
    if col_start == ncols - 1 
        % Assume that the last row containes the secular varation coefficients
        dcoeff = table2array(mag_tab(:, col_start + 1)); 
    else 
        % Need to manually compute the 5-years coefficient variation
        next_coeff = table2array(mag_tab(:, col_start + 1)); 
        dcoeff = 0.2*(next_coeff - coeff);
    end

    % Maximum available model degree
    maxdeg = max(deg);

    % Note: the order goes from (0, deg) so it has one more coefficient wrt the degree.

    % Sine (h) and cosine (g) coefficients
    g = zeros(maxdeg, maxdeg + 1);
    h = zeros(maxdeg, maxdeg + 1);

    % Secular variations of the coefficients 
    dg = zeros(maxdeg, maxdeg + 1);
    dh = zeros(maxdeg, maxdeg + 1);

    % Fill the coefficient matrices
    for i = 1:length(letter)

        degi = deg(i); 
        ordi = ord(i) + 1; % 1-Based notation 

        if strcmp(letter{i}, 'g') 
            g(degi, ordi) = coeff(i); 
            dg(degi, ordi) = dcoeff(i); 
        else 
            h(degi, ordi) = coeff(i); 
            dh(degi, ordi) = dcoeff(i);
        end

    end

    % Matrix to store the transformation factors from the Schmidt to the
    % Gauss normalisation. The Gauss functions can then be computed as 
    % P(Gauss) = S * P(Schmidt). These factors are efficient because they 
    % are independent on the point coordinates and thus must be computed just once.  
    S = zeros(maxdeg, maxdeg + 1);

    % Matrices to store the Gaussian normalised IGRF coefficients!
    % The columns are: [degree, order, gn, dgn, hn, dhn]
    ncoeffs = maxdeg*(maxdeg+1)/2 + maxdeg;
    out_mat = zeros(ncoeffs, 6); 
    
    idx = 1;
    for n = 1:maxdeg 
        for m = 0:n 

            IDm = m + 1; % 1-Based index

            if m >= 2   
                Snm = S(n, m)*sqrt((n-m+1)/(n+m));
            elseif m == 1
                Snm = S(n, m)*sqrt(2*n/(n+1));
            elseif n > 1 % n > 1 and m == 0
                Snm = S(n-1, 1)*(2*n-1)/n;
            else % n == 1 and m == 0 
                Snm = 1; 
            end

            S(n, IDm) = Snm;    
            
            % Fills the matrix
            out_mat(idx, 1) = n; 
            out_mat(idx, 2) = m; 
            out_mat(idx, 3) = Snm*g(n, IDm);
            out_mat(idx, 4) = Snm*dg(n, IDm);
            out_mat(idx, 5) = Snm*h(n, IDm);
            out_mat(idx, 6) = Snm*dh(n, IDm);

            idx = idx + 1;
        end
    end

    % Write the two matrix to the desired text file.
    dlmwrite(sprintf("%s.txt", output_name), out_mat, '\t');

end