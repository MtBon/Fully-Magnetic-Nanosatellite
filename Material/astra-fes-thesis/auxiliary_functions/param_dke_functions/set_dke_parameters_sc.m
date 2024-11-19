function set_dke_parameters_sc(sim_type)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Set Parameters for DKE of ASTRA FES
    %
    % (c) ASTRA - Aerospace Science and Technology Dept. - PoliMi -
    % 2018
    %
    % Requires the simulation type: OS or MC to set dispersion
    % Others needed parameters can be recovered from the calling function as
    % PARAM=evalin('caller','PARAM');
    %
    % V 2.0
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Parameters Set Up
    if nargin == 0
        sim_type = 'OS';
    end

    % Set dispersions
    switch sim_type
        case 'OS'
            sc.errI = 0.1/3;        % 10% - 3σ
            sc.errCOM = 0;
            sc.errdif = 0.25/3;     % 25% - 3σ
            sc.errspe = 0.25/3;     % 25% - 3σ
            sc.errdrag = 0.1/3;     % 10% - 3σ
            
        case 'MC'
            sc.errI = 0.25/3;       % 25% - 3σ
            sc.errCOM = 5/1000;     % 5mm - Uniform
            sc.errdif = 0.5/3;      % 50% - 3σ
            sc.errspe = 0.5/3;      % 50% - 3σ
            sc.errdrag = 0.25/3;    % 25% - 3σ
        otherwise
            error('ERROR: sim type must be ''OS'' or ''MC'' ')
    end

    %% Flags to Define Panels Type
    sc.panels_open_flag = 1;    % 1 = wing panels are open 
    sc.bol_flag = 0;            % 0 = end of life - 1 = begin of life

    %% Spacecraft Properties
    
    % Mass (kg)
    sc.m_sc = 16.5;
    
    % Mean S/C Area (m²):
    sc.A = ((0.1*0.3) + ((0.4 + 0.1*sqrt(2))*0.3))/2; 

    % Inverse of the Ballistic Coefficient (m²/kg)
    sc.BC = single(2.2*sc.A/6.42); 

    % Inertia Matrix with the wheels (kg*m²)
    sc.Ixx =  0.066709;
    sc.Iyy =  0.066894;
    sc.Izz =  0.025517;
    sc.Ixy = -0.008119;
    sc.Ixz = -0.000064;
    sc.Iyz = -0.000419;
    
    % Error on Inertia Matrix
    sc.rIxx = (1 + sc.errI*randn(1))*sc.Ixx;
    sc.rIyy = (1 + sc.errI*randn(1))*sc.Iyy;
    sc.rIzz = (1 + sc.errI*randn(1))*sc.Izz;
    sc.rIxy = (1 + sc.errI*randn(1))*sc.Ixy;
    sc.rIxz = (1 + sc.errI*randn(1))*sc.Ixz;
    sc.rIyz = (1 + sc.errI*randn(1))*sc.Iyz;
    
    % Matrix Assembly
    sc.I = [sc.rIxx, sc.rIxy, sc.rIxz;
            sc.rIxy, sc.rIyy, sc.rIyz; 
            sc.rIxz, sc.rIyz, sc.rIzz];
         
    sc.invI = inv(sc.I); %#ok<*NASGU>

    % C.o.M. Offset with respect to Geometric Reference Frame (IB).
    sc.Xcom = 50.000/1000; % (m)
    sc.Ycom = 48.070/1000; % (m)
    sc.Zcom = 177.88/1000; % (m)
    
    % Error on C.O.M - Uniform Distribution
    sc.COM = [sc.Xcom + 2*sc.errCOM*(rand(1) - 0.5); 
              sc.Ycom + 2*sc.errCOM*(rand(1) - 0.5); 
              sc.Zcom + 2*sc.errCOM*(rand(1) - 0.5)];
           

    %% Data for Surfaces.
    % Number of Faces
    sc.NS = 10;
    
    % Geometric dimensions of the body
    L_s = 0.1; % Small side (m)
    L_l = 0.3; % Long side (m)
    L_p = 0.2; % Panel Extension (m)

    % Normalized Normal Directions of the Faces - In the IB frame
    N1_1 = [ 1;  0;  0]; 
    N1_2 = [-1;  0;  0]; 
    N1_3 = [ 0;  1;  0]; 
    N1_4 = [ 0; -1;  0]; 
    N1_5 = [ 0;  0;  1]; 
    N1_6 = [ 0;  0; -1]; 
    
    % Solar Panels Normal Directions - In the IB Frame 
    N1_7  = [ -1; -1; 0]/sqrt(2); % Positive X (front SA)
    N1_8  = [  1;  1; 0]/sqrt(2); % Positive X
    N1_9  = [ -1; -1; 0]/sqrt(2); % Positive Y (front SA)
    N1_10 = [  1;  1; 0]/sqrt(2); % Positive Y

    sc.N1 = [N1_1, N1_2, N1_3, N1_4, N1_5, N1_6, N1_7, N1_8, N1_9, N1_10];
    
    % Center of pressure vectors of the Faces with respect to CoM in the IB Frame.
    R1_1 = [   L_s; L_s/2; L_l/2] - sc.COM;
    R1_2 = [  0.00; L_s/2; L_l/2] - sc.COM;
    R1_3 = [ L_s/2;   L_s; L_l/2] - sc.COM;
    R1_4 = [ L_s/2;  0.00; L_l/2] - sc.COM;
    R1_5 = [ L_s/2; L_s/2;   L_l] - sc.COM;
    R1_6 = [ L_s/2; L_s/2;  0.00] - sc.COM;
    
    % Center of pressure vectors of the Solar Panels in the IB Frame.
    R1_7   = [ L_s + L_p/2/sqrt(2);     - L_p/2/sqrt(2); L_l/2] - sc.COM;
    R1_8   = [ L_s + L_p/2/sqrt(2);     - L_p/2/sqrt(2); L_l/2] - sc.COM;
    R1_9   = [     - L_p/2/sqrt(2); L_s + L_p/2/sqrt(2); L_l/2] - sc.COM;
    R1_10  = [     - L_p/2/sqrt(2); L_s + L_p/2/sqrt(2); L_l/2] - sc.COM;
    
    sc.R1 = [R1_1, R1_2, R1_3, R1_4, R1_5, R1_6, R1_7, R1_8, R1_9, R1_10];

    % Faces Area subdivision
    S1_sv = [0, 0, 0, 0, 1, 1, 0, 0, 0, 0]; % Small Faces
    S1_lv = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0]; % Large Faces
    S1_pv = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1]; % Panels Faces
    
    % Faces Area (m²)
    S1_s = L_s*L_s;
    S1_l = L_s*L_l;
    S1_p = L_l*L_p;
    
    sc.S1 = S1_s.*S1_sv+S1_l.*S1_lv+S1_p.*S1_pv;

    sc.R_surf = sc.R1;
    sc.N_surf = sc.N1;
    sc.S_surf = sc.S1;

    %% If accurate Solar shadowing model is needed (SRP/Solar Arrays) Uncomment This
    % NOTE: you sould uncomment also the related block in the Simulink SRP function
    %
    % Surfaces triangulation for self shadowing
    % [TRI,face_id,is_panel,tri_num_face,tri_num_panel,XYZ,TRI_center,TRI_normal] = wings_triangulation(NS,L_s,L_l,L_p);

    %% Data for Solar Panels
    sc.P0_Sun = 1367; % W/m2
    sc.S_panel = 0.002651;  %m2 - Single cell area
    
    % BOL vs EOL
    if sc.bol_flag
        sc.eta_panel = 0.296; % Panels Efficiency (BOL)
        sc.Id = 0.988179; % Degradation Factor - Determines 1.06W per cell (BOL)
    else
        sc.eta_panel = 0.279; % Panels Efficiency (EOL-5e14)
        sc.Id = 0.935639; % Degradation Factor - Determines 0.946W per cell (EOL-5e14)
    end
    
    % Solar Cells Distribution
    if sc.panels_open_flag
        sc.S_panels = sc.S_panel*[2,4,0,4,0,2,14,0,14,0]; % m2 - Number of cells per face
    else
        sc.S_panels = sc.S_panel*[2,7,0,7,0,2,0,0,0,0]; % #ok<*UNRCH> %m2 - Number of cells per face
    end
    
    %% Data for Solar pressure
    sc.c_dif = 0.1;
    sc.c_spe = 0.1;

    % Coefficients computation 
    sc.c_dif = (1 + sc.errdif*randn(1))*sc.c_dif;
    sc.c_spe = (1 + sc.errspe*randn(1))*sc.c_spe;
    sc.c_abs = 1 - sc.c_dif - sc.c_spe;


    %% Atmospheric Drag Data
    % C_D for a flat plate perpendicular to flow (2D flow)
    sc.C_D = 2.2;
    sc.C_D = (1 + sc.errdrag*randn(1))*sc.C_D;

    %% Parasitic Magnetic Dipole Data
    % Magnetic Dipole
    d = randn(3, 1); % Am2
    d = d./norm(d);
    sc.d = (sc.m_sc*2.5e-3).*d;

    % Eddy Current coefficient : Adapted From NASA (1969a).
    % It uses a thin spherical shell of radius (re), thickness (te) and conductivity (se) 
    sc.te = 0.005; % (m)
    sc.re = 0.05;  % (m)
    sc.se = 37.7e6; % Al conductivity [S/m]
    sc.Keddy = 2/3*pi*sc.re^4*sc.se*sc.te;

    %% Data for Flexibility - Only used in FLEX DKE simulink model
    if evalin('caller','flex_flag')
        % Flex modelled as a small sloshing in all the S/C
        sc.flex.rhoXe=5.894; %Xe Density kg/m3
        sc.flex.Dt=0.1; %Tank Diameter m
        sc.flex.Ht=0.3; %Tank Height m
        sc.flex.Vthank=pi*sc.flex.Dt*sc.flex.Dt*sc.flex.Ht/4; %Thank Volume m3
        sc.flex.ml=sc.flex.rhoXe*sc.flex.Vthank; %kg
        % Slosh Parameters from slosh modelling
        sc.flex.m1=(sc.flex.ml*(sc.flex.Dt/2)/(2.2*sc.flex.Ht))*tanh(1.84*sc.flex.Ht/(sc.flex.Dt/2)); %kg - sloshing mass
        sc.flex.h1=sc.flex.Ht/2-((sc.flex.Dt/2)/(1.84))*tanh(1.84*sc.flex.Ht/(sc.flex.Dt/2)); %m - sloshing position
        sc.flex.m0=sc.flex.ml-sc.flex.m1; %kg - fix sloshing mass
        sc.flex.h0=(sc.flex.ml/sc.flex.m0)*(sc.flex.Ht/2-((sc.flex.Dt/2)^2/(2*sc.flex.Ht)))-sc.flex.h1*sc.flex.m1/sc.flex.m0; %m - fix sloshing position
        sc.flex.Kc=sc.flex.m1*(9.81/1.19/sc.flex.Ht)*((tanh(1.84*sc.flex.Ht/(sc.flex.Dt/2)))^2); %[kg m2/s2] - Equivalent spring coefficient
        sc.flex.Ke=sc.flex.Kc; %[kg m2/s2] - Equivalent spring coefficient
        sc.flex.z1=0.05; %Damping Ratio
        sc.flex.Cc=2*sc.flex.m1*sc.flex.z1*sqrt(sc.flex.Kc/sc.flex.m1); %Equivalent damping coefficient [kg m2/s]
        sc.flex.Ce=sc.flex.Cc; %Equivalent damping coefficient [kg m2/s]
    end

    %% Export SC variable to the caller workspace
    assignin('caller', 'sc', eval('sc'))

return
