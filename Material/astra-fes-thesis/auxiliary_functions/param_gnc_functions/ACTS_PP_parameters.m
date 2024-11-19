function [TUN_STRUCT,NTUN_STRUCT,SAMPLE_TIME] = ACTS_PP_parameters
% it returns a struct array containing tunable struct parameters
% for ACTS_PP Block Models

% Sample time
SAMPLE_TIME.time_1=1; %[s]
% ... add other sample times if present

% Tunable
TUN_STRUCT.var_1=single(1);
TUN_STRUCT.var_2=double(2);
TUN_STRUCT.var_3=uint32(3);
TUN_STRUCT.var_4=single(69);


% Static - Non Tunable
NTUN_STRUCT.var_1=single(33);
NTUN_STRUCT.var_2=boolean(22);
NTUN_STRUCT.var_3=single(48);
NTUN_STRUCT.var_4=single(69);

% ....

end
