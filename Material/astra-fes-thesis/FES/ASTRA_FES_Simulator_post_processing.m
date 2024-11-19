clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics - Orbit Attitude with GNC and Environment Simulator -
% Post Processing Results - V2.0
%
%
%
% © MAT - Aerospace Science and Technology Dept. - PoliMi - 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../AuxiliaryFiles'))
addpath(genpath('../EnvironmentData'))
addpath(genpath('../../DataBase'))

%% Load File
[filename_out, path_out] = uigetfile('*.mat', 'Please select an output data file');
load(fullfile(path_out,filename_out))


%% ----------- Postprocessing Single Shot analisys -----------
if strcmp(sim_type,'OS')
    
    %% Initialize Variables
    n_step_sol = max(size(output{1}.t_vect));
    n_step_sol_adcs = max(size(output{1}.t_adcs));
    
    loading_vector = round(linspace(1,n_step_sol,n_step_sol_adcs)); % vector containing the position of the considered points
    loading_vector_adcs = 1:1:n_step_sol_adcs; % adcs vector containing the position of the considered points
    
    %% Extract Output
    t_vect=output{1}.t_vect(loading_vector);
    t_adcs = output{1}.t_adcs(loading_vector_adcs);
    r=output{1}.r(loading_vector,:);
    v=output{1}.v(loading_vector,:);
    q=output{1}.q(loading_vector,:);
    w=output{1}.w(loading_vector,:);
    ww=output{1}.ww(loading_vector,:);
    f_pert=output{1}.f_pert(loading_vector,:);
    m_pert=output{1}.m_pert(loading_vector,:);
    q_ref=output{1}.q_ref(loading_vector_adcs,:);
    Sc2S=output{1}.Sc2S(loading_vector,:);
    Ecl=output{1}.Ecl(loading_vector);
    McW=output{1}.McW(loading_vector_adcs,:);
    try
        McW_actuated=output{1}.McW_actuated(loading_vector,:);
    catch
        McW_actuated=output{1}.McW_actuated(loading_vector_adcs,:);
    end
    DcMag=output{1}.DcMag(loading_vector_adcs,:);
    try
        DcMag_actuated=output{1}.DcMag_actuated(loading_vector,:);
    catch
        DcMag_actuated=output{1}.DcMag_actuated(loading_vector_adcs,:);
    end
    PanelsPower=output{1}.PanelsPower(loading_vector,:);
    PanelsTotalPower=output{1}.PanelsTotalPower(loading_vector);
    PanelsTotalEnergy=output{1}.PanelsTotalEnergy(loading_vector);
    
    %% Compute OS Derived Values
    x_error = zeros(1,n_step_sol_adcs);
    y_error = zeros(1,n_step_sol_adcs);
    z_error = zeros(1,n_step_sol_adcs);
    x_axis = zeros(3,n_step_sol_adcs);
    y_axis = zeros(3,n_step_sol_adcs);
    z_axis = zeros(3,n_step_sol_adcs);
    
    roll = zeros(1,n_step_sol_adcs);
    pitch= zeros(1,n_step_sol_adcs);
    yaw  = zeros(1,n_step_sol_adcs);
    
    % Compute DCM Matrices and Euler angles
    for jj = 1:n_step_sol_adcs
        % DCM Matrix
        CI_= [1-2*q(jj,2)^2-2*q(jj,3)^2, 2*(q(jj,1)*q(jj,2)-q(jj,3)*q(jj,4)), 2*(q(jj,3)*q(jj,1)+q(jj,2)*q(jj,4));...
            2*(q(jj,1)*q(jj,2)+q(jj,3)*q(jj,4)), 1-2*q(jj,3)^2-2*q(jj,1)^2, 2*(q(jj,2)*q(jj,3)-q(jj,1)*q(jj,4));...
            2*(q(jj,3)*q(jj,1)-q(jj,2)*q(jj,4)), 2*(q(jj,2)*q(jj,3)+q(jj,1)*q(jj,4)), 1-2*q(jj,1)^2-2*q(jj,2)^2];
        
        % Body axes
        x_axis(:,jj) =  CI_(:,1);
        y_axis(:,jj) =  CI_(:,2);
        z_axis(:,jj) =  CI_(:,3);
        
        % Reference DCM Matrix
        CI_ref= [1-2*q_ref(jj,2)^2-2*q_ref(jj,3)^2, 2*(q_ref(jj,1)*q_ref(jj,2)-q_ref(jj,3)*q_ref(jj,4)), 2*(q_ref(jj,3)*q_ref(jj,1)+q_ref(jj,2)*q_ref(jj,4));...
            2*(q_ref(jj,1)*q_ref(jj,2)+q_ref(jj,3)*q_ref(jj,4)), 1-2*q_ref(jj,3)^2-2*q_ref(jj,1)^2, 2*(q_ref(jj,2)*q_ref(jj,3)-q_ref(jj,1)*q_ref(jj,4));...
            2*(q_ref(jj,3)*q_ref(jj,1)-q_ref(jj,2)*q_ref(jj,4)), 2*(q_ref(jj,2)*q_ref(jj,3)+q_ref(jj,1)*q_ref(jj,4)), 1-2*q_ref(jj,1)^2-2*q_ref(jj,2)^2];
        
        % Pointing error
        x_error(jj) = acosd(dot(x_axis(:,jj),CI_ref(:,1)));
        y_error(jj) = acosd(dot(y_axis(:,jj),CI_ref(:,2)));
        z_error(jj) = acosd(dot(z_axis(:,jj),CI_ref(:,3)));
        
        % LVLH frame
        rr = r(jj,:)./norm(r(jj,:));
        vv = v(jj,:)./norm(v(jj,:));
        R_bar = rr;
        H_bar = cross(rr,vv);
        V_bar = cross(H_bar,R_bar);
        
        ROT_lvlh = [R_bar; V_bar; H_bar];
        A_bl = CI_*ROT_lvlh';
        x_lvlh = A_bl(:,1);
        y_lvlh = A_bl(:,2);
        z_lvlh = A_bl(:,3);
        
        z_ort_h = z_lvlh - (dot(z_lvlh,H_bar)/norm(H_bar)^2)*H_bar';
        z_ort_v = z_lvlh - (dot(z_lvlh,V_bar)/norm(V_bar)^2)*V_bar';
        x_ort_v = x_lvlh - (dot(x_lvlh,V_bar)/norm(V_bar)^2)*V_bar';
        y_ort_r = y_lvlh - (dot(y_lvlh,R_bar)/norm(R_bar)^2)*R_bar';
        
        roll(jj) = acosd(dot(R_bar,z_ort_v));
        %roll(jj) = acosd(dot(-H_bar,x_ort_v));
        pitch(jj)= acosd(dot(R_bar,z_ort_h));
        yaw(jj)  = acosd(dot(V_bar,y_ort_r));
    end
    
    %% Compute Power Consumption Values
    % Magnetorquer values
    Iz=3.3/25.5; %I=V/R at 0.34Am2 as Datasheet
    Ixy=3.3/28.5; %I=V/R at 0.31Am2 as Datasheet
    nsz=0.34/Iz; %Proportional coefficient of magnetorquer: Dc=nsI
    nsxy=0.31/Ixy; %Proportional coefficient of magnetorquer: Dc=nsI
    PMag=[28.5*(DcMag_actuated(:,1)./nsxy).^2, 28.5*(DcMag_actuated(:,2)./nsxy).^2, 25.5*(DcMag_actuated(:,3)./nsz).^2]; %W - Magnetorquer Power
    PMagTotal=sum(PMag,2);
    
    % Wheels values
    q0=interp1([0,0.25e-3,0.5e-3,1e-3,1.5e-3],[0,0.04,0.07,0.12,0.28],abs(McW_actuated)); %A at w0
    m0=interp1([0,0.25e-3,0.5e-3,1e-3,1.5e-3],[1.27e-4,2.86e-4,3.48e-4,4.88e-4,4.88e-4],abs(McW_actuated)); %A/(rad/s)
    I=q0+abs(ww).*m0; %A
    PWheels=1.1*5*I; %W - 5V with 10% margin
    PWheelsTotal=sum(PWheels,2);
    
    %% Plotting Script
    OS_plotting;
    
    %% ----------- Postprocessing MonteCarlo analisys -----------
elseif strcmp(sim_type,'MC')
    %% Prompt fraction of point to plot
    prompt = 'Insert the point fraction to be analyzed ( from 0 to 1 ; 1 for all points): ';
    point_fraction = input(prompt);
    
    %% Initialize MC Variables
    n_step_sol_tot = max(size(output{1}.q));            % number of total points
    n_step_sol = round(n_step_sol_tot*point_fraction);  % fraction of points considered in the pp analisys
    n_step_sol_adcs = max(size(output{1}.t_adcs));      % number of total adcs points
    n_step_pp = min(n_step_sol,n_step_sol_adcs);        % number of points considered for the post processing analysis
    
    %% 'com' number of points control
    if strcmp(adcs_mode_type,'com')
        if n_step_pp < n_step_sol_adcs
            n_step_pp = n_step_sol_adcs;
            disp('Fraction of point too low \n')
            fprintf('Fraction of point considered: %f \n', n_step_sol_adcs/n_step_sol_tot)
        end
    end
    
    loading_vector = round(linspace(1,n_step_sol_tot,n_step_pp)); % vector containing the position of the considered points
    loading_vector_adcs = 1:1:n_step_sol_adcs; % adcs vector containing the position of the considered points
    
    Q0 = zeros(4,n_sim);   % initial quaternions
    W0 = zeros(3,n_sim);   % initial angular velocity
    WW0 = zeros(4,n_sim);  % initial wheels velocities
    t_vect = zeros(n_step_pp,n_sim); % time vector
    t_adcs = zeros(n_step_sol_adcs,n_sim); % time vector
    r  = zeros(n_step_pp,3,n_sim);   % ECI radius [m]
    q  = zeros(n_step_pp,4,n_sim);   % quaternions
    w  = zeros(n_step_pp,3,n_sim);   % angular velocity
    ww = zeros(n_step_pp,4,n_sim);   % wheel ang velocities
    f_pert = zeros(n_step_pp,3,n_sim);     % perturbation force
    m_pert = zeros(n_step_pp,3,n_sim);     % perturbations torque
    q_ref =  zeros(n_step_sol_adcs,4,n_sim);     % reference quaternion
    Sc2S =   zeros(n_step_pp,3,n_sim);     % spacecraft to sun vector
    Sc2GS = zeros(n_step_sol_adcs-1,3,n_sim); % spacecraft to ground station
    Ecl = zeros(n_step_pp,n_sim);          % Eclipse flag
    McW = zeros(n_step_sol_adcs,4,n_sim);   % Wheels control torque
    McW_actuated = zeros(n_step_pp,4,n_sim);   %Actuated Wheels control torque
    DcMag = zeros(n_step_sol_adcs,3,n_sim); % Magnetotorquers dipole
    DcMag_actuated = zeros(n_step_pp,3,n_sim); %Actuated Magnetotorquers dipole
    ThermalSun = zeros(n_step_pp,NS,n_sim);    % Thermal power from sun
    ThermalEarth = zeros(n_step_pp,NS,n_sim);  % Thermal power from Earth
    ThermalPower = zeros(n_step_pp,NS,n_sim);  % Total thermal power
    x_axis = zeros(3,n_step_pp); % x axis dir in ECI
    y_axis = zeros(3,n_step_pp);  % y axis dir in ECI
    z_axis = zeros(3,n_step_pp);  % z axis dir in ECI
    x_error = zeros(n_step_pp,n_sim); % x axis pointing error
    y_error = zeros(n_step_pp,n_sim); % y axis pointing error
    z_error = zeros(n_step_pp,n_sim); % z axis pointing error
    PanelsTotalPower = zeros(n_step_pp,1,n_sim); % Panels power from sun
    PMag=zeros(n_step_pp,3,n_sim); %Magnetorquer Power
    PMagTotal=zeros(n_step_pp,1,n_sim); %Magnetorquer Total Power
    PWheels=zeros(n_step_pp,4,n_sim); %Wheels Power
    PWheelsTotal=zeros(n_step_pp,1,n_sim); %Wheels Total Power
    PTotal=zeros(n_step_pp,1,n_sim);
    
    %% Extract MC Output
    for ii = 1:n_sim
        % Loading data from output variable
        Q0(:,ii)=output{ii}.Q0;
        W0(:,ii)=output{ii}.W0;
        WW0(:,ii)=output{ii}.WW0;
        t_vect(:,ii)=output{ii}.t_vect(loading_vector);
        t_adcs(:,ii)=output{ii}.t_adcs(loading_vector_adcs);
        r(:,:,ii)=output{ii}.r(loading_vector,:);
        q(:,:,ii)=output{ii}.q(loading_vector,:);
        w(:,:,ii)=output{ii}.w(loading_vector,:);
        ww(:,:,ii)=output{ii}.ww(loading_vector,:);
        f_pert(:,:,ii)=output{ii}.f_pert(loading_vector,:);
        m_pert(:,:,ii)=output{ii}.m_pert(loading_vector,:);
        q_ref(:,:,ii)=output{ii}.q_ref(loading_vector_adcs,:);
        Sc2S(:,:,ii)=output{ii}.Sc2S(loading_vector,:);
        if strcmp(adcs_mode_type,'com')
            Sc2GS(:,:,ii)=output{ii}.Sc2GS(loading_vector_adcs(1:n_step_sol_adcs-1),:);
        end
        Ecl(:,ii)=output{ii}.Ecl(loading_vector,:);
        McW(:,:,ii)=output{ii}.McW(loading_vector_adcs,:);
        McW_actuated(:,:,ii)=output{ii}.McW_actuated(loading_vector,:);
        DcMag(:,:,ii)=output{ii}.DcMag(loading_vector_adcs,:);
        DcMag_actuated(:,:,ii)=output{ii}.DcMag_actuated(loading_vector,:);
        ThermalSun(:,:,ii) = output{ii}.ThermalSun(loading_vector,:);
        ThermalEarth(:,:,ii) = output{ii}.ThermalEarthSimple(loading_vector,:);
        ThermalPower(:,:,ii) = output{ii}.ThermalPower(loading_vector,:);
        PanelsTotalPower(:,:,ii)=output{ii}.PanelsTotalPower(loading_vector);
        
        %% Compute MC Derived Value
        % Compute DCM Matrices
        for jj = 1:n_step_pp
            % Computing axis directions from q
            x_axis(:,jj) = [ 1-2*q(jj,2,ii)^2-2*q(jj,3,ii)^2 ; 2*(q(jj,1,ii)*q(jj,2,ii)+q(jj,3,ii)*q(jj,4,ii)) ; 2*(q(jj,3,ii)*q(jj,1,ii)-q(jj,2,ii)*q(jj,4,ii)) ];
            y_axis(:,jj) = [ 2*(q(jj,1,ii)*q(jj,2,ii)-q(jj,3,ii)*q(jj,4,ii)) ; 1-2*q(jj,3,ii)^2-2*q(jj,1,ii)^2 ;  2*(q(jj,2,ii)*q(jj,3,ii)+q(jj,1,ii)*q(jj,4,ii)) ];
            z_axis(:,jj) = [  2*(q(jj,3,ii)*q(jj,1,ii)+q(jj,2,ii)*q(jj,4,ii)) ; 2*(q(jj,2,ii)*q(jj,3,ii)-q(jj,1,ii)*q(jj,4,ii)) ; 1-2*q(jj,1,ii)^2-2*q(jj,2,ii)^2 ];
            
            % Reference DCM Matrix
            CI_ref= [1-2*q_ref(jj,2,ii)^2-2*q_ref(jj,3,ii)^2, 2*(q_ref(jj,1,ii)*q_ref(jj,2,ii)-q_ref(jj,3,ii)*q_ref(jj,4,ii)), 2*(q_ref(jj,3,ii)*q_ref(jj,1,ii)+q_ref(jj,2,ii)*q_ref(jj,4,ii));...
                2*(q_ref(jj,1,ii)*q_ref(jj,2,ii)+q_ref(jj,3,ii)*q_ref(jj,4,ii)), 1-2*q_ref(jj,3,ii)^2-2*q_ref(jj,1,ii)^2, 2*(q_ref(jj,2,ii)*q_ref(jj,3,ii)-q_ref(jj,1,ii)*q_ref(jj,4,ii));...
                2*(q_ref(jj,3,ii)*q_ref(jj,1,ii)-q_ref(jj,2,ii)*q_ref(jj,4,ii)), 2*(q_ref(jj,2,ii)*q_ref(jj,3,ii)+q_ref(jj,1,ii)*q_ref(jj,4,ii)), 1-2*q_ref(jj,1,ii)^2-2*q_ref(jj,2,ii)^2];
            
            % Computing pointing errors
            x_error(jj,ii) = acosd(x_axis(1,jj)*CI_ref(1,1) + x_axis(2,jj)*CI_ref(2,1) + x_axis(3,jj)*CI_ref(3,1));
            y_error(jj,ii) = acosd(y_axis(1,jj)*CI_ref(1,2) + y_axis(2,jj)*CI_ref(2,2) + y_axis(3,jj)*CI_ref(3,2));
            z_error(jj,ii) = acosd(z_axis(1,jj)*CI_ref(1,3) + z_axis(2,jj)*CI_ref(2,3) + z_axis(3,jj)*CI_ref(3,3));
        end
        
        % Compute Power Consumption Values
        % Magnetorquer values
        Iz=3.3/25.5; %I=V/R at 0.34Am2 as Datasheet
        Ixy=3.3/28.5; %I=V/R at 0.31Am2 as Datasheet
        nsz=0.34/Iz; %Proportional coefficient of magnetorquer: Dc=nsI
        nsxy=0.31/Ixy; %Proportional coefficient of magnetorquer: Dc=nsI
        PMag(:,:,ii)=[28.5*(DcMag_actuated(:,1,ii)./nsxy).^2, 28.5*(DcMag_actuated(:,2,ii)./nsxy).^2, 25.5*(DcMag_actuated(:,3,ii)./nsz).^2]; %W - Magnetorquer Power
        PMagTotal(:,1,ii)=sum(PMag(:,:,ii),2);
        
        % Wheels values
        q0=interp1([0,0.25e-3,0.5e-3,1e-3,1.5e-3],[0,0.04,0.07,0.12,0.28],abs(McW_actuated(:,:,ii))); %A at w0
        m0=interp1([0,0.25e-3,0.5e-3,1e-3,1.5e-3],[1.27e-4,2.86e-4,3.48e-4,4.88e-4,4.88e-4],abs(McW_actuated(:,:,ii))); %A/(rad/s)
        I=q0+abs(ww(:,:,ii)).*m0; %A
        PWheels(:,:,ii)=1.1*5*I; %W - 5V with 10% margin
        PWheelsTotal(:,1,ii)=sum(PWheels(:,:,ii),2);
        
        % ADCS total power
        PTotal(:,1,ii)=PWheelsTotal(:,1,ii)+PMagTotal(:,1,ii);
    end
    clear output
    %% Switch Case Plot
    switch adcs_mode_type
        
        case 'det'
            
            MC_DET_plotting;
            
        case 'des'
            
            MC_DES_plotting;
            
        case 'poi'
            
            MC_POINT_plotting;
            
        case 'sle'
            
            MC_SLE_plotting;
            
        case 'sun'
            
            MC_SUN_plotting;
            
        case 'com'
            
            MC_COM_plotting;
            
        case 'zen'
            
            MC_ZEN_plotting;
            
        case 'seq'
            
            %To be implemented!
            
        otherwise
            error('ADCS mode not valid');
    end
    t = t_vect(:,1);
end
% cleanfigure('minimumPointsDistance', 0.1)
% matlab2tikz('test.tex',...
%             'width', '5.5cm',...
%             'extraAxisOptions',['font=\tiny,'...
%                                 'xlabel style={font=\scriptsize},'...
%                                 'ylabel style={font=\scriptsize},'...
%                                 'zlabel style={font=\scriptsize},'...
%                                 'legend style={font=\scriptsize},'])
%myFigTemplate
