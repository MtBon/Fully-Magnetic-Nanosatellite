function GPS=HP_GPS_init(sensor_time)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPS errors initialization %  
% GOM SPACE - GPS KITS Sensor 
% NovAtel OEM719 GPS module
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Time [s]
% Sample time output dello strumento
GPS.sampletime = sensor_time;              %[s] sample time
% Delay (150% 3-sigma) - GPS Delay any output
GPS.delay = 0.015*(ones(1,1)+1.5/3*randn(1,1)); % 15 ms
GPS.delay_noise = 0.001/3*(ones(1,1)+1.5/3*randn(1,1)); % 1 ms (3-sigma)
% Noise (10% 3-sigma) 
GPS.noiseT = 20e-9*(ones(1,1)+0.1/3*randn(1,1));  % 20 ns  m (1-sigma)
GPS.noiseseedT = floor(2928845*rand(1,1)); % Numero grande a caso

%% Position [m]
% Bias (10% 3-sigma) 
GPS.bias = 0.1*(sign(randn(1,3)).*ones(1,3)+0.1/3*randn(1,3)); % 10cm - Within the inbox of the spacecraft
% Noise (10% 3-sigma) 
GPS.noise = 1.5*(ones(1,3)+0.1/3*randn(1,3));  % 1.5 m (1-sigma)
GPS.noiseseed = floor(4925845*rand(1,3)); % Numero grande a caso
% Risoluzione (gain standard)
GPS.resolution = 1.5/2^12;           % 12-Bit LSB per range +/-1.5 m

%% Velocity [m/s]
% Bias
GPS.biasV = zeros(1,3);   % 0 m
% Noise (10% 3-sigma) 
GPS.noiseV = 0.03*(ones(1,3)+0.1/3*randn(1,3));    % 0.03 m/s (1-sigma)
GPS.noiseseedV = floor(293675845*rand(1,3)); % Numero grande a caso
% Velocity Time Latency (50% 3-sigma) [s]
GPS.V_latency=0.1*(ones(1,1)+0.5/3*randn(1,1)); % 0.1 s 
GPS.V_latency_noise=0.005/3*(ones(1,1)+1.5/3*randn(1,1)); % 5 ms (3-sigma)
% Risoluzione (gain standard)
GPS.resolutionV = 0.03/2^12;      % 12-Bit LSB per range +/-0.03 m/s

