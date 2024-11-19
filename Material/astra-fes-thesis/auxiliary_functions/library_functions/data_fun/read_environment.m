function env = read_environment(filename_env, sim_start, sim_duration)

    % Retrieve the environment parameters to be used for the simulation
    % settings and dynamical propagations. 
    %
    % Parameters
    % ----------
    %   filename_env: string
    %       Name of the file including the simulation time settings.
    %   sim_start: double
    %       Julian Date at the start of the simulation
    %   sim_duration: double
    %       Simulation duration, in seconds.
    %
    % Returns 
    % -------
    %   env: struct
    %       Environment data structure

    % Load settings from file
    data = read_json(filename_env);

    % Starting Julian Date and Time settings
    env.jdate0 = sim_start;
    env.jdate0_int = floor(sim_start);
    env.jdate0_frac = sim_start - env.jdate0_int;

    % Load IERS EOP Data 
    env.eop = parse_eop(data.EOP.file, sim_start, sim_start + sim_duration/86400);

    % Retrieve leapseconds count at the simulation start time
    LEAP = load_leapseconds(data.Leapseconds.file); 
    env.leap = get_leapseconds(sim_start - 2451545, LEAP);

    % Load Magnetic Model Data
    env.mag_data = load_mag_model(data.MagneticModel.file); 
    
    % Load Gravity Model Data
    env.grav_data = load_gravity_model(data.GravityModel.file, data.GravityModel.degree);

    % Gravitational Parameters in (m³/s²)
    env.GM_EARTH = env.grav_data.GM;
    env.GM_MOON = 4.902800076000000e+12;
    env.GM_SUN = 1.327124400409440e+20;

    % Earth and Sun radii (m)
    env.RADIUS_EARTH = 1e3*6378.137;
    env.RADIUS_SUN = 1e3*696000;
    
    % Earth Atmosphere Altitude (m). This value is used to 
    % account for the atmosphere thickness when computing the 
    % shadow coefficient.
    env.ATMO_HEIGHT = 90000;

    % Earth J2 Coefficient (-)
    env.J2 = 0.00108263;
    
    % Earth Air Density model look-up table (km, kg/m³)
    env.rho_h = density_exponential_model();

    % View Factor Table Coefficients
    env.vf.table = reshape([...
        0.9946, 0.9752, 0.9062, 0.8282, 0.6962, 0.5922, 0.4432, 0.2490, 0.1110, ... 
        0.0279, 0.9594, 0.9303, 0.8523, 0.7755, 0.6521, 0.5560, 0.4177, 0.2346, ...
        0.1040, 0.0260, 0.8636, 0.8174, 0.7215, 0.6438, 0.5332, 0.4526, 0.3398, ...
        0.1917, 0.0854, 0.0214, 0.7986, 0.7460, 0.6428, 0.5638, 0.4575, 0.3840, ...
        0.2853, 0.1604, 0.0718, 0.0181, 0.7246, 0.6676, 0.5600, 0.4809, 0.3792, ...
        0.3122, 0.2267, 0.1248, 0.0554, 0.0139, 0.6438, 0.5839, 0.4752, 0.3983, ...
        0.3030, 0.2427, 0.1692, 0.0881, 0.0375, 0.0092, 0.5585, 0.4973, 0.3906, ...
        0.3182, 0.2318, 0.1790, 0.1173, 0.0547, 0.0209, 0.0047, 0.4713, 0.4107, ...
        0.3090, 0.2430, 0.1673, 0.1230, 0.0737, 0.0285, 0.0085, 0.0015, 0.3852, ...
        0.3271, 0.2331, 0.1747, 0.1112, 0.0762, 0.0401, 0.0113, 0.0016, 0.0000, ...
        0.3025, 0.2487, 0.1650, 0.1156, 0.0654, 0.0403, 0.0172, 0.0026, 0.0000, ...
        0.0000, 0.2256, 0.1774, 0.1065, 0.0677, 0.0320, 0.0163, 0.0045, 0.0000, ...
        0.0000, 0.0000, 0.1570, 0.1158, 0.0597, 0.0326, 0.0112, 0.0036, 0.0000, ...
        0.0000, 0.0000, 0.0000, 0.0991, 0.0665, 0.0267, 0.0107, 0.0012, 0.0000, ...
        0.0000, 0.0000, 0.0000, 0.0000, 0.0216, 0.0091, 0, 0, 0, 0, 0, 0, 0, 0, ...
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 10, 15);

    env.vf.bp1 = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 5];
    env.vf.bp2 = [0, 20, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 160, 180];
    
end