function IMU=HP_IMU_MS3025_init(sensor_time)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMU errors initialization % 
% MEMSENSE 3025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gyros operating on deg/s units
% Accelerometers operating on G units

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note on Noise Amplitude and Noise Density
% Noise Density is measured in U/sqrt(Hz) --> (A/sqrt(samplerate))^2=D
% Noise Amplitude in U --> sqrt(D*samplerate)=A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variance Limits
IMU.gyros.max_var=1e-1;
IMU.accelerometers.max_var=1e-5;

% Sampling Frequencies
IMU.gyros.samplerate = 1/sensor_time;
IMU.accelerometers.samplerate = 1/sensor_time;

% Measurement Range
% 75 deg/s
% +/-4 g 
IMU.gyros.range = 75*ones(1,3);
IMU.accelerometers.range = 4*ones(1,3);

% Resolution
% 32 bits - single precision floating (2^23)
IMU.gyros.resolution = IMU.gyros.range/2^23;
IMU.accelerometers.resolution = IMU.accelerometers.range/2^23;

% Static Bias
% 10.8 deg/h (10% 3-sigma)
% 270 ug (600% 3-sigma - max 2000ug)
IMU.gyros.bias = (10.8/3600)*(sign(randn(1,3)).*ones(1,3)+0.1/3*randn(1,3));
IMU.accelerometers.bias = 270e-6*(sign(randn(1,3)).*ones(1,3)+6/3*randn(1,3));

% White Noise (Stochastic Process) - Angle random walk for gyroscopes or Velocity random walk for accelerometers
% D=0.0038 deg/S/sqrt(Hz) (30% 3-sigma)
% D=20.4 ug/sqrt(Hz) (30% 3-sigma)
% Tuning Parameters Included
IMU.gyros.noise = 0.75*0.0038*(ones(1,3)+0.3/3*randn(1,3));
IMU.accelerometers.noise = 0.75*20.4e-6*(ones(1,3)+0.3/3*randn(1,3));

% Random Walk (Stochastic Process) - Rate random walk for gyroscopes or Acceleration random walk
% D=0.15 deg/h (30% 3-sigma)
% D=0.008 m/s/h (30% 3-sigma)
% Tuning Parameters Included
IMU.gyros.randomwalk = 0.5*0.15/3600*(ones(1,3)+0.3/3*randn(1,3));
IMU.accelerometers.randomwalk = 0.1*0.008/3600*(ones(1,3)+0.3/3*randn(1,3));

% Bias Stability (Stochastic Process)
% 0.96 deg/h (10% 3-sigma -correlation time 100s)
% 3.7 ug (10% 3-sigma - correlation time 100s)
IMU.gyros.BS = 0.96/3600*(ones(1,3)+0.1/3*randn(1,3));
IMU.accelerometers.BS = 3.7e-6*(ones(1,3)+0.1/3*randn(1,3));
IMU.gyros.ctime = sensor_time*ones(1,3);
IMU.accelerometers.ctime = sensor_time*ones(1,3);

% Bias Stability and Random Walk Seeds
IMU.gyros.noiseseed = ceil(1000*rand(1,3));
IMU.gyros.BSseed = ceil(1000*rand(1,3));
IMU.accelerometers.noiseseed = ceil(1000*rand(1,3));
IMU.accelerometers.BSseed = ceil(1000*rand(1,3));

% Scale Factor Errors (Stochastic)
% SFx=x*SFx*1e-6;
% Gyros Various (ppm) (100% 3-sigma)
% Accelerometers Various (ppm) (400% 3-sigma)
IMU.gyros.SFE = [1263, 1263, 1087].*(ones(1,3)+1/3*randn(1,3));
IMU.accelerometers.SFE = [34.34, 34.34, 185.65].*(ones(1,3)+4/3*randn(1,3));
% Gyros 10 ppm (10% 3-sigma)
% Accelerometers 1 ppm (10% 3-sigma)
IMU.gyros.SFA = 10*(ones(1,3)+0.1/3*randn(1,3));
IMU.accelerometers.SFA = 1*(ones(1,3)+0.1/3*randn(1,3));
% Gyros 500 ppm (0.05%) (10% 3-sigma)
% Accelerometers 3000 ppm (0.3%) (10% 3-sigma)
IMU.gyros.SFN = (500*(ones(1,3)+0.1/3*randn(1,3)))./IMU.gyros.range(1);
IMU.accelerometers.SFN = (3000*(ones(1,3)+0.1/3*randn(1,3)))./IMU.accelerometers.range(1);

% Axes Cross-sensitivity - Not symmetric
% 2% (3-sigma)
IMU.gyros.mismat = eye(3) + 0.02/3*[0 randn randn; randn 0 randn; randn randn 0];
IMU.accelerometers.mismat = eye(3) + 0.02/3*[0 randn randn; randn 0 randn; randn randn 0];

% IMU Location w.r.t. CoG in Body Frame - [m]
IMU.gyros.PI_location = [0 0 -0.1];
IMU.accelerometers.PI_location = IMU.gyros.PI_location;

% Body to IMU Axes Rotation Matrix      
IMU.gyros.mountmat = eye(3); % Transform PI to IMU frame (IMU_A_PI)
IMU.accelerometers.mountmat = IMU.gyros.mountmat;

% IMU Error Rotation Matrix 
% Random Axis and 0.1deg (3sigma)
ax_err=randn(1,3);
IMU.gyros.mounterr = vrrotvec2mat([ax_err./norm(ax_err),deg2rad(0.1/3*randn)]);
IMU.accelerometers.mounterr = IMU.gyros.mounterr;


