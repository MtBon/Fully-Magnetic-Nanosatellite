function output_compressor(workspace, frac_points, compress_adcs)

    load(workspace); %#ok<LOAD>

    % Total number of simulation points
    n_step_sol_tot = length(output{1}.t_vect);  %#ok<*NODEF> 
    
    % Total number of points considered in the post-process analysis
    n_step_sol = round(n_step_sol_tot*frac_points); 
    
    % Indexes of the considered points
    idx = round(linspace(1, n_step_sol_tot, n_step_sol));
    
    % Total number of points for the ADCS
    n_step_sol_tot_adcs = length(output{1}.t_adcs);  %#ok<*NODEF> 

    if compress_adcs == 1
        % Indexes of the considered ADCS variables
        idx_adcs = round(linspace(1, n_step_sol_tot_adcs, n_step_sol));
    else
        idx_adcs = 1:1:n_step_sol_tot_adcs;
    end

    output2 = cell(n_sim,1);
    
    for i = 1:n_sim
        
        % DKE Variables --------------------
        output2{i}.Q0                  = output{i}.Q0;
        output2{i}.W0                  = output{i}.W0;
        output2{i}.WW0                 = output{i}.WW0;
        output2{i}.oevi                = output{i}.oevi;
        output2{i}.t_vect              = output{i}.t_vect(idx);
        output2{i}.r                   = output{i}.r(idx, :);
        output2{i}.v                   = output{i}.v(idx, :);
        output2{i}.a                   = output{i}.a(idx, :);
        output2{i}.q                   = output{i}.q(idx, :);
        output2{i}.q_ref               = output{i}.q_ref(idx, :);
        output2{i}.Sc2S                = output{i}.Sc2S(idx, :);
        output2{i}.Ecl                 = output{i}.Ecl(idx, :);
        output2{i}.Alb                 = output{i}.Alb(idx, :);
        output2{i}.w                   = output{i}.w(idx, :);
        output2{i}.ww                  = output{i}.ww(idx, :);
        output2{i}.b                   = output{i}.b(idx, :);
        output2{i}.f_pert              = output{i}.f_pert(idx, :);
        output2{i}.m_pert              = output{i}.m_pert(idx, :);
        
        % ADCS variables -----
        output2{i}.t_adcs              = output{i}.t_adcs(idx_adcs);
        output2{i}.McW                 = output{i}.McW(idx_adcs, :);
        output2{i}.DcMag               = output{i}.DcMag(idx_adcs, :);
        output2{i}.McW_actuated        = output{i}.McW_actuated(idx_adcs, :);
        output2{i}.DcMag_actuated      = output{i}.DcMag_actuated(idx_adcs, :);
        
        output2{i}.McW_actuated        = output{i}.McW_actuated(idx, :);
        output2{i}.DcMag_actuated      = output{i}.DcMag_actuated(idx, :);
        output2{i}.PanelsPower         = output{i}.PanelsPower(idx, :);
        output2{i}.PanelsTotalPower    = output{i}.PanelsTotalPower(idx, :);
        output2{i}.PanelsTotalEnergy   = output{i}.PanelsTotalEnergy(idx, :);
        output2{i}.ThermalPower        = output{i}.ThermalPower(idx, :);
        output2{i}.ThermalSun          = output{i}.ThermalSun(idx, :);
        output2{i}.ThermalEarthSimple  = output{i}.ThermalEarthSimple (idx, :);

        output{i} = output2{i}; %#ok<AGROW>
    end

    clear output2;
    new_workspace = strcat(workspace,'_compressed');
    save(new_workspace,'-v7.3');
    
end

