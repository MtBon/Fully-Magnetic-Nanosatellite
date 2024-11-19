%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Montecarlo ('MC') ZENITH Mode Simulation Plotting Script
% Script of plots called by DKE_ACS_MAIN_Simulator_post_processing
%
% DO NOT RUN ALONE!
%
% © MAT - Aerospace Science and Technology Dept. - PoliMi - 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Restrict computation to first 6 body surfaces
% N.B. MUST be modified to make these computations with NS>6
NS=6;

%% COM Computation
% Initialize variables
earth_angle =zeros(n_step_pp,NS,n_sim); % Earth Angle
x_axis = zeros(3,n_step_pp); % x axis dir in ECI
y_axis = zeros(3,n_step_pp);  % y axis dir in ECI
z_axis = zeros(3,n_step_pp);  % z axis dir in ECI

for ii = 1:n_sim
    % t_sim loop
    for jj = 1:n_step_pp
        % Computing earth (nadir) angles
        x_axis(:,jj) = [ 1-2*q(jj,2,ii)^2-2*q(jj,3,ii)^2 ; 2*(q(jj,1,ii)*q(jj,2,ii)+q(jj,3,ii)*q(jj,4,ii)) ; 2*(q(jj,3,ii)*q(jj,1,ii)-q(jj,2,ii)*q(jj,4,ii)) ];
        y_axis(:,jj) = [ 2*(q(jj,1,ii)*q(jj,2,ii)-q(jj,3,ii)*q(jj,4,ii)) ; 1-2*q(jj,3,ii)^2-2*q(jj,1,ii)^2 ;  2*(q(jj,2,ii)*q(jj,3,ii)+q(jj,1,ii)*q(jj,4,ii)) ];
        z_axis(:,jj) = [  2*(q(jj,3,ii)*q(jj,1,ii)+q(jj,2,ii)*q(jj,4,ii)) ; 2*(q(jj,2,ii)*q(jj,3,ii)-q(jj,1,ii)*q(jj,4,ii)) ; 1-2*q(jj,1,ii)^2-2*q(jj,2,ii)^2 ];
            
        earth_angle(jj,1,ii) = acosd( -( x_axis(1,jj)*r(jj,1,ii) + x_axis(2,jj)*r(jj,2,ii) + x_axis(3,jj)*r(jj,3,ii)) /norm(r(jj,:,ii)));
        earth_angle(jj,2,ii) = acosd( -(-x_axis(1,jj)*r(jj,1,ii) - x_axis(2,jj)*r(jj,2,ii) - x_axis(3,jj)*r(jj,3,ii)) /norm(r(jj,:,ii)));
        earth_angle(jj,3,ii) = acosd( -( y_axis(1,jj)*r(jj,1,ii) + y_axis(2,jj)*r(jj,2,ii) + y_axis(3,jj)*r(jj,3,ii)) /norm(r(jj,:,ii)));
        earth_angle(jj,4,ii) = acosd( -(-y_axis(1,jj)*r(jj,1,ii) - y_axis(2,jj)*r(jj,2,ii) - y_axis(3,jj)*r(jj,3,ii)) /norm(r(jj,:,ii)));
        earth_angle(jj,5,ii) = acosd( -( z_axis(1,jj)*r(jj,1,ii) + z_axis(2,jj)*r(jj,2,ii) + z_axis(3,jj)*r(jj,3,ii)) /norm(r(jj,:,ii)));
        earth_angle(jj,6,ii) = acosd( -(-z_axis(1,jj)*r(jj,1,ii) - z_axis(2,jj)*r(jj,2,ii) - z_axis(3,jj)*r(jj,3,ii)) /norm(r(jj,:,ii)));
    end
end

%% ZEN Plots
%Error Trend
figure(111)
plot(t_vect(:,1)/3600,x_error(:,1),'r',t_vect(:,1)/3600,y_error(:,1),'g',t_vect(:,1)/3600,z_error(:,1),'b')
grid on
hold on
for ii = 2:n_sim
    plot(t_vect(:,ii)/3600,x_error(:,ii),'r',t_vect(:,ii)/3600,y_error(:,ii),'g',t_vect(:,ii)/3600,z_error(:,ii),'b')
end
plot(t_vect(:,1)/3600,3+zeros(size(x_error)),'-.k')
xlabel('Time [h]')
ylabel('Error [deg]')
legend('x-axis','y-axis','z-axis','3 deg. Requirement')

%Earth Angle Error Trend
figure(1)
plot(t_vect(1:end,1)/3600,earth_angle(1:end,1,1),'r',t_vect(1:end,1)/3600,earth_angle(1:end,2,1),'g',t_vect(1:end,1)/3600,earth_angle(1:end,3,1),'b',t_vect(1:end,1)/3600,earth_angle(1:end,4,1),'y',t_vect(1:end,1)/3600,earth_angle(1:end,5,1),'c',t_vect(1:end,1)/3600,earth_angle(1:end,6,1),'m')
grid on
hold on
for ii = 2:n_sim
    plot(t_vect(1:end,ii)/3600,earth_angle(1:end,1,ii),'r',t_vect(1:end,ii)/3600,earth_angle(1:end,2,ii),'g',t_vect(1:end,ii)/3600,earth_angle(1:end,3,ii),'b',t_vect(1:end,ii)/3600,earth_angle(1:end,4,ii),'y',t_vect(1:end,ii)/3600,earth_angle(1:end,5,ii),'c',t_vect(1:end,ii)/3600,earth_angle(1:end,6,ii),'m')
end
hold on
xlabel('Time [h]')
ylabel('Earth Angles [deg]')
legend('+x','-x','+y','-y','+z','-z')

%Wheels Angular Velocity
figure(2)
plot(t_vect(1:end,1)/3600,ww(1:end,1,1),'r',t_vect(1:end,1)/3600,ww(1:end,2,1),'g',t_vect(1:end,1)/3600,ww(1:end,3,1),'b',t_vect(1:end,1)/3600,ww(1:end,4,1),'y')
grid on
hold on
plot(t_vect(1:end,1)/3600,0*ww(1:end,1,1)-0.1*RWL.MaxWheelW,'-.r','linewidth',1.5)
plot(t_vect(1:end,1)/3600,0*ww(1:end,1,1)-0.7*RWL.MaxWheelW,'-.k','linewidth',1.5)
plot(t_vect(1:end,1)/3600,0*ww(1:end,1,1)+0.7*RWL.MaxWheelW,'-.k','linewidth',1.5)
plot(t_vect(1:end,1)/3600,0*ww(1:end,1,1)+0.1*RWL.MaxWheelW,'-.r','linewidth',1.5)
for ii = 2:n_sim
    plot(t_vect(1:end,ii)/3600,ww(1:end,1,ii),'r',t_vect(1:end,ii)/3600,ww(1:end,2,ii),'g',t_vect(1:end,ii)/3600,ww(1:end,3,ii),'b',t_vect(1:end,ii)/3600,ww(1:end,4,ii),'y')
end
xlabel('t [h]')
ylabel('Wheel Angular Velocity [rad/s]')
legend('\omega_{A}', '\omega_{B}', '\omega_{C}','\omega_{D}','10 % Saturation','70 % Saturation')

%Angular Velocity
figure(3)
plot(t_vect(1:end,1)/3600,w(1:end,1,1),'r',t_vect(1:end,1)/3600,w(1:end,2,1),'g',t_vect(1:end,1)/3600,w(1:end,3,1),'b')
hold on
grid on
plot(t_vect(1:end,1)/3600, (w(1:end,1).^2 + w(1:end,2).^2 + w(1:end,3).^2).^0.5,'-.c' )
plot(t_vect(1:end,1)/3600, 0*w(1:end,3) + w_treshold,'--k')
plot(t_vect(1:end,1)/3600, 0*w(1:end,3) - w_treshold,'--k')
for ii = 2:n_sim
    plot(t_vect(1:end,ii)/3600,w(1:end,1,ii),'r',t_vect(1:end,ii)/3600,w(1:end,2,ii),'g',t_vect(1:end,ii)/3600,w(1:end,3,ii),'b')
    plot(t_vect(1:end,ii)/3600, (w(1:end,1,ii).^2 + w(1:end,2,ii).^2 + w(1:end,3,ii).^2).^0.5,'-.c' )
end
xlabel('t [h]')
ylabel('Angular Velocity [rad/s]')
legend('\omega_x', '\omega_y', '\omega_z','Magnitude','Treshold')

%Control Torque (Wheel)
figure(4)
plot(t_adcs(1:end,1)/3600,McW(1:end,1,1),'r',t_adcs(1:end,1)/3600,McW(1:end,2,1),'g',t_adcs(1:end,1)/3600,McW(1:end,3,1),'b',t_adcs(1:end,1)/3600,McW(1:end,4,1),'y')
grid on
hold on
plot(t_adcs(1:end,1)/3600,0*McW(1:end,1,1)-1*RWL.MaxTorque,'-.r','linewidth',1.5)
plot(t_adcs(1:end,1)/3600,0*McW(1:end,1,1)+1*RWL.MaxTorque,'-.r','linewidth',1.5)
for ii = 2:n_sim
    plot(t_adcs(1:end,ii)/3600,McW(1:end,1,ii),'r',t_adcs(1:end,ii)/3600,McW(1:end,2,ii),'g',t_adcs(1:end,ii)/3600,McW(1:end,3,ii),'b',t_adcs(1:end,ii)/3600,McW(1:end,4,ii),'y')
end
xlabel('t [h]')
ylabel('Wheel Control Torque [N]')
legend('Mw_1', 'Mw_2', 'Mw_3', 'Mw_4','Max Torque')

%Actuated Control Torque (Wheel)
figure(5)
plot(t_vect(1:end,1)/3600,McW_actuated(1:end,1,1),'r',t_vect(1:end,1)/3600,McW_actuated(1:end,2,1),'g',t_vect(1:end,1)/3600,McW_actuated(1:end,3,1),'b',t_vect(1:end,1)/3600,McW_actuated(1:end,4,1),'y')
grid on
hold on
plot(t_adcs(1:end,1)/3600,0*McW(1:end,1,1)-1*RWL.MaxTorque,'-.r','linewidth',1.5)
plot(t_adcs(1:end,1)/3600,0*McW(1:end,1,1)+1*RWL.MaxTorque,'-.r','linewidth',1.5)
for ii = 2:n_sim
    plot(t_vect(1:end,ii)/3600,McW_actuated(1:end,1,ii),'r',t_vect(1:end,ii)/3600,McW_actuated(1:end,2,ii),'g',t_vect(1:end,ii)/3600,McW_actuated(1:end,3,ii),'b',t_vect(1:end,ii)/3600,McW_actuated(1:end,4,ii),'y')
end
xlabel('t [h]')
ylabel('Actuated Wheel Control Torque [N]')
legend('Mw_1', 'Mw_2', 'Mw_3', 'Mw_4','Max Torque')

%Wheels power consumption
figure(6)
plot(t_vect(:,1)/3600,PWheelsTotal(:,1,1),'k')
hold on
for ii = 2:n_sim
    plot(t_vect(:,ii)/3600,PWheelsTotal(:,1,ii),'k')
end
xlabel('Time  [h]')
ylabel('Wheel Power Consumption [W]')
legend('P_{w_{tot}}')
grid on

%ADCS actuators power consumption
figure(7)
plot(t_vect(:,1)/3600,PTotal(:,1,1),'k')
hold on
for ii = 2:n_sim
    plot(t_vect(:,ii)/3600,PTotal(:,1,ii),'k')
end
xlabel('Time  [h]')
ylabel('ADCS actuator power consumption [W]')
legend('P_{ADCS_{tot}}')
grid on

%Panels Total Power
figure(8)
plot(t_vect(:,1)/3600,PanelsTotalPower(:,1,1),'k',t_vect(:,1)/3600,mean(PanelsTotalPower(:,1,1))*ones(size(t_vect(:,1))),'r')
hold on
plot(t_vect(:,1)/3600,mean(PanelsTotalPower(PanelsTotalPower(:,1,1)~=0,1,1))*ones(size(PanelsTotalPower(:,1,1))),'Color',[0.5,0.8,0])
for ii = 2:n_sim
    plot(t_vect(:,ii)/3600,PanelsTotalPower(:,1,ii),'k',t_vect(:,ii)/3600,mean(PanelsTotalPower(:,1,ii))*ones(size(t_vect(:,ii))),'r')
    plot(t_vect(:,ii)/3600,mean(PanelsTotalPower(PanelsTotalPower(:,1,ii)~=0,1,ii))*ones(size(PanelsTotalPower(:,1,ii))),'Color',[0.5,0.8,0])
end
xlabel('Time [h]')
ylabel('Generated Power [W]')
legend('P_{tot}','P_{avg}','P_{avg_{Sun}}')
grid on

%% FigTemplate
FigTemplateMC