%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Single Shot ('OS') Simulation Plotting Script
% Script of plots called by DKE_ACS_MAIN_Simulator_post_processing
%
% DO NOT RUN ALONE!
%
% © MAT - Aerospace Science and Technology Dept. - PoliMi - 2019
%
% V 2.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General Plots
%Quaternions
figure(1)
plot(t_vect(1:end)/3600,q(1:end,1),t_vect(1:end)/3600,q(1:end,2),t_vect(1:end)/3600,q(1:end,3),t_vect(1:end)/3600,q(1:end,4))
xlabel('t [h]')
ylabel('Quaternions[-]')
legend('q_{v_1}', 'q_{v_2}', 'q_{v_3}','q_s');
grid on
hold on

figure(21)
plot(t_vect/3600,roll,'r','Linewidth',1)
hold on
grid on
plot(t_vect/3600,pitch,'b','Linewidth',1)
plot(t_vect/3600,yaw,'k','Linewidth',1)
xlabel('Time [h]')
ylabel('Axis angle [deg]')
legend('roll','pitch','yaw')

%Angular Velocity
figure(2)
plot(t_vect(1:end)/3600,w(1:end,1),t_vect(1:end)/3600,w(1:end,2),t_vect(1:end)/3600,w(1:end,3))
hold on
grid on
plot(t_vect(1:end)/3600, (w(1:end,1).^2 + w(1:end,2).^2 + w(1:end,3).^2).^0.5,'-.' )
plot(t_vect(1:end)/3600, 0*w(1:end,3) + w_treshold,'--k')
plot(t_vect(1:end)/3600, 0*w(1:end,3) - w_treshold,'--k')
xlabel('t [h]')
ylabel('Angular Velocity [rad/s]')
legend('\omega_x', '\omega_y', '\omega_z','Magnitude','Treshold')

%Orbit plot
figure(3)
plot3(r(:,1)./1000,r(:,2)./1000,r(:,3)./1000,'r')
hold on
grid on
axis equal
set(gca, 'NextPlot','add','Color','none');
[xx, yy, zz] = ellipsoid(0, 0, 0, req, req, req, 72);
earth = surf(xx, yy, -zz);
cdata = imread('earth.jpg');
set(earth, 'FaceColor', 'texturemap', 'CData', cdata,'FaceAlpha', 0.75, 'EdgeColor', 'none');
colormap('gray');
xlabel('x [km]', 'FontSize', 12);
ylabel('y [km]', 'FontSize', 12);
zlabel('z [km]', 'FontSize', 12);

%Wheels Angular Velocity
figure(4)
plot(t_vect(1:end)/3600,ww(1:end,1),t_vect(1:end)/3600,ww(1:end,2),t_vect(1:end)/3600,ww(1:end,3),t_vect(1:end)/3600,ww(1:end,4))
grid on
hold on
plot(t_vect(1:end)/3600,0*ww(1:end,1)-0.1*RWL.MaxWheelW,'-.r','linewidth',1.5)
plot(t_vect(1:end)/3600,0*ww(1:end,1)-0.7*RWL.MaxWheelW,'-.k','linewidth',1.5)
plot(t_vect(1:end)/3600,0*ww(1:end,1)+0.7*RWL.MaxWheelW,'-.k','linewidth',1.5)
plot(t_vect(1:end)/3600,0*ww(1:end,1)+0.1*RWL.MaxWheelW,'-.r','linewidth',1.5)
xlabel('t [h]')
ylabel('Wheel Angular Velocity [rad/s]')
legend('\omega_{A}', '\omega_{B}', '\omega_{C}','\omega_{D}','10 % Saturation','70 % Saturation')

%F perturbations
figure(5)
plot(t_vect(1:end)/3600,f_pert(1:end,1),t_vect(1:end)/3600,f_pert(1:end,2),t_vect(1:end)/3600,f_pert(1:end,3))
xlabel('t [h]')
ylabel('Perturbation Force [N]')
legend('f_x', 'f_y', 'f_z')
grid on
hold on

%M perturbations
figure(6)
plot(t_vect(1:end)/3600,m_pert(1:end,1),t_vect(1:end)/3600,m_pert(1:end,2),t_vect(1:end)/3600,m_pert(1:end,3))
xlabel('t [h]')
ylabel('Perturbation Torque [Nm]')
legend('m_x', 'm_y', 'm_z')
grid on
hold on

%Control Torque (Wheel)
figure(7)
plot(t_adcs/3600,McW(1:end,1),t_adcs/3600,McW(1:end,2),t_adcs/3600,McW(1:end,3),t_adcs/3600,McW(1:end,4))
xlabel('t [h]')
ylabel('Wheel Control Torque [N]')
legend('Mw_1', 'Mw_2', 'Mw_3', 'Mw_4')
grid on
hold on

%Actuated Control Torque (Wheel)
figure(8)
plot(t_vect/3600,McW_actuated(1:end,1),t_vect/3600,McW_actuated(1:end,2),t_vect/3600,McW_actuated(1:end,3),t_vect/3600,McW_actuated(1:end,4))
xlabel('t [h]')
ylabel('Actuated Wheel Control Torque [N]')
legend('Mw_1', 'Mw_2', 'Mw_3', 'Mw_4')
grid on
hold on


%Control Dipole (MagTorq)
figure(9)
plot(t_adcs/3600,DcMag(1:end,1),t_adcs/3600,DcMag(1:end,2),t_adcs/3600,DcMag(1:end,3))
xlabel('t [h]')
ylabel('Magnetorquers Dipole [Am2]')
legend('Dm_x', 'Dm_y', 'Dm_z')
grid on
hold on

%Actuated Control Dipole (MagTorq)
figure(10)
plot(t_vect/3600,DcMag_actuated(1:end,1),t_vect/3600,DcMag_actuated(1:end,2),t_vect/3600,DcMag_actuated(1:end,3))
xlabel('t [h]')
ylabel('Actuated Magnetorquers Dipole [Am2]')
legend('Dm_x', 'Dm_y', 'Dm_z')
grid on
hold on

%% ADCS Mode Based Plots
%Pointing Error
if strcmp(adcs_mode_type,'poi') || strcmp(adcs_mode_type,'sle') || strcmp(adcs_mode_type,'com') || strcmp(adcs_mode_type,'zen')
    figure(19)
    plot(t_vect/3600,x_error,'r','Linewidth',1)
    hold on
    grid on
    plot(t_vect/3600,y_error,'b','Linewidth',1)
    plot(t_vect/3600,z_error,'k','Linewidth',1)
    plot(t_vect/3600,3+x_error*0,'-.b','Linewidth',1.5)
    xlabel('Time [h]')
    ylabel('Error [deg]')
    legend('x-axis','y-axis','z-axis','3 deg. Requirement')
end

%Sun Angles
if strcmp(adcs_mode_type,'sun') || strcmp(adcs_mode_type,'sle')
    % Initialize variables
    sun_angle =zeros(n_step_sol_adcs,NS);
    for jj = 1:n_step_sol_adcs
        % Computing sun angles
        sun_angle(jj,1) = acosd( ( x_axis(1,jj)*Sc2S(jj,1) + x_axis(2,jj)*Sc2S(jj,2) + x_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        sun_angle(jj,2) = acosd( (-x_axis(1,jj)*Sc2S(jj,1) - x_axis(2,jj)*Sc2S(jj,2) - x_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        sun_angle(jj,3) = acosd( ( y_axis(1,jj)*Sc2S(jj,1) + y_axis(2,jj)*Sc2S(jj,2) + y_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        sun_angle(jj,4) = acosd( (-y_axis(1,jj)*Sc2S(jj,1) - y_axis(2,jj)*Sc2S(jj,2) - y_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        sun_angle(jj,5) = acosd( ( z_axis(1,jj)*Sc2S(jj,1) + z_axis(2,jj)*Sc2S(jj,2) + z_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        sun_angle(jj,6) = acosd( (-z_axis(1,jj)*Sc2S(jj,1) - z_axis(2,jj)*Sc2S(jj,2) - z_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
    end
    figure(11)
    plot(t_vect(:,1)/3600,sun_angle(:,1,1),t_vect(:,1)/3600,sun_angle(:,2,1),t_vect(:,1)/3600,sun_angle(:,3,1),t_vect(:,1)/3600,sun_angle(:,4,1),t_vect(:,1)/3600,sun_angle(:,5,1),t_vect(:,1)/3600,sun_angle(:,6,1))
    hold on
    grid on
    xlabel('Time [h]')
    ylabel('Sun Angles [deg]')
    legend('+x','-x','+y','-y','+z','-z')
end


%Sun and Earth Angles
if strcmp(adcs_mode_type,'com')
    % Initialize variables
    Sc2GS=output{1}.Sc2GS(loading_vector_adcs(1:n_step_sol_adcs-1),:);
    sun_angle =zeros(n_step_sol_adcs,NS);
    gs_angle =zeros(n_step_sol_adcs,NS);
    earth_angle =zeros(n_step_sol_adcs,NS);
    for jj = 1:n_step_sol_adcs
        % Computing sun angles
        sun_angle(jj,1) = acosd( ( x_axis(1,jj)*Sc2S(jj,1) + x_axis(2,jj)*Sc2S(jj,2) + x_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        sun_angle(jj,2) = acosd( (-x_axis(1,jj)*Sc2S(jj,1) - x_axis(2,jj)*Sc2S(jj,2) - x_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        sun_angle(jj,3) = acosd( ( y_axis(1,jj)*Sc2S(jj,1) + y_axis(2,jj)*Sc2S(jj,2) + y_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        sun_angle(jj,4) = acosd( (-y_axis(1,jj)*Sc2S(jj,1) - y_axis(2,jj)*Sc2S(jj,2) - y_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        sun_angle(jj,5) = acosd( ( z_axis(1,jj)*Sc2S(jj,1) + z_axis(2,jj)*Sc2S(jj,2) + z_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        sun_angle(jj,6) = acosd( (-z_axis(1,jj)*Sc2S(jj,1) - z_axis(2,jj)*Sc2S(jj,2) - z_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        
        % Computing earth (nadir) angles
        earth_angle(jj,1) = acosd( -( x_axis(1,jj)*r(jj,1) + x_axis(2,jj)*r(jj,2) + x_axis(3,jj)*r(jj,3)) /norm(r(jj,:)));
        earth_angle(jj,2) = acosd( -(-x_axis(1,jj)*r(jj,1) - x_axis(2,jj)*r(jj,2) - x_axis(3,jj)*r(jj,3)) /norm(r(jj,:)));
        earth_angle(jj,3) = acosd( -( y_axis(1,jj)*r(jj,1) + y_axis(2,jj)*r(jj,2) + y_axis(3,jj)*r(jj,3)) /norm(r(jj,:)));
        earth_angle(jj,4) = acosd( -(-y_axis(1,jj)*r(jj,1) - y_axis(2,jj)*r(jj,2) - y_axis(3,jj)*r(jj,3)) /norm(r(jj,:)));
        earth_angle(jj,5) = acosd( -( z_axis(1,jj)*r(jj,1) + z_axis(2,jj)*r(jj,2) + z_axis(3,jj)*r(jj,3)) /norm(r(jj,:)));
        earth_angle(jj,6) = acosd( -(-z_axis(1,jj)*r(jj,1) - z_axis(2,jj)*r(jj,2) - z_axis(3,jj)*r(jj,3)) /norm(r(jj,:)));
    end
        
    for jj = 1:n_step_sol_adcs-1   
        % Computing Antenna angles
        gs_angle(jj,1) = acosd( ( x_axis(1,jj)*Sc2GS(jj,1) + x_axis(2,jj)*Sc2GS(jj,2) + x_axis(3,jj)*Sc2GS(jj,3)) /norm(Sc2GS(jj,:)));
        gs_angle(jj,2) = acosd( (-x_axis(1,jj)*Sc2GS(jj,1) - x_axis(2,jj)*Sc2GS(jj,2) - x_axis(3,jj)*Sc2GS(jj,3)) /norm(Sc2GS(jj,:)));
        gs_angle(jj,3) = acosd( ( y_axis(1,jj)*Sc2GS(jj,1) + y_axis(2,jj)*Sc2GS(jj,2) + y_axis(3,jj)*Sc2GS(jj,3)) /norm(Sc2GS(jj,:)));
        gs_angle(jj,4) = acosd( (-y_axis(1,jj)*Sc2GS(jj,1) - y_axis(2,jj)*Sc2GS(jj,2) - y_axis(3,jj)*Sc2GS(jj,3)) /norm(Sc2GS(jj,:)));
        gs_angle(jj,5) = acosd( ( z_axis(1,jj)*Sc2GS(jj,1) + z_axis(2,jj)*Sc2GS(jj,2) + z_axis(3,jj)*Sc2GS(jj,3)) /norm(Sc2GS(jj,:)));
        gs_angle(jj,6) = acosd( (-z_axis(1,jj)*Sc2GS(jj,1) - z_axis(2,jj)*Sc2GS(jj,2) - z_axis(3,jj)*Sc2GS(jj,3)) /norm(Sc2GS(jj,:)));        
    end
    
    figure(16)
    plot(t_vect(:,1)/3600,sun_angle(:,1,1),t_vect(:,1)/3600,sun_angle(:,2,1),t_vect(:,1)/3600,sun_angle(:,3,1),t_vect(:,1)/3600,sun_angle(:,4,1),t_vect(:,1)/3600,sun_angle(:,5,1),t_vect(:,1)/3600,sun_angle(:,6,1))
    hold on
    grid on
    xlabel('Time [h]')
    ylabel('Sun Angles [deg]')
    legend('+x','-x','+y','-y','+z','-z')
    
    figure(17)
    plot(t_adcs(:,1)/3600,gs_angle(:,1,1),t_adcs(:,1)/3600,gs_angle(:,2,1),t_adcs(:,1)/3600,gs_angle(:,3,1),t_adcs(:,1)/3600,gs_angle(:,4,1),t_adcs(:,1)/3600,gs_angle(:,5,1),t_adcs(:,1)/3600,gs_angle(:,6,1))
    hold on
    grid on
    xlabel('Time [h]')
    ylabel('GS Angles [deg]')
    legend('+x','-x','+y','-y','+z','-z')
    
    figure(18)
    plot(t_vect(:,1)/3600,earth_angle(:,1,1),t_vect(:,1)/3600,earth_angle(:,2,1),t_vect(:,1)/3600,earth_angle(:,3,1),t_vect(:,1)/3600,earth_angle(:,4,1),t_vect(:,1)/3600,earth_angle(:,5,1),t_vect(:,1)/3600,earth_angle(:,6,1))
    hold on
    grid on
    xlabel('Time [h]')
    ylabel('Earth Angles [deg]')
    legend('+x','-x','+y','-y','+z','-z')
end

if strcmp(adcs_mode_type,'zen')
    % Initialize variables
    earth_angle =zeros(n_step_sol_adcs,NS);
    for jj = 1:n_step_sol_adcs
        % Computing sun angles
        sun_angle(jj,1) = acosd( ( x_axis(1,jj)*Sc2S(jj,1) + x_axis(2,jj)*Sc2S(jj,2) + x_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        sun_angle(jj,2) = acosd( (-x_axis(1,jj)*Sc2S(jj,1) - x_axis(2,jj)*Sc2S(jj,2) - x_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        sun_angle(jj,3) = acosd( ( y_axis(1,jj)*Sc2S(jj,1) + y_axis(2,jj)*Sc2S(jj,2) + y_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        sun_angle(jj,4) = acosd( (-y_axis(1,jj)*Sc2S(jj,1) - y_axis(2,jj)*Sc2S(jj,2) - y_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        sun_angle(jj,5) = acosd( ( z_axis(1,jj)*Sc2S(jj,1) + z_axis(2,jj)*Sc2S(jj,2) + z_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        sun_angle(jj,6) = acosd( (-z_axis(1,jj)*Sc2S(jj,1) - z_axis(2,jj)*Sc2S(jj,2) - z_axis(3,jj)*Sc2S(jj,3)) /norm(Sc2S(jj,:)));
        
        % Computing earth (nadir) angles
        earth_angle(jj,1) = acosd( -( x_axis(1,jj)*r(jj,1) + x_axis(2,jj)*r(jj,2) + x_axis(3,jj)*r(jj,3)) /norm(r(jj,:)));
        earth_angle(jj,2) = acosd( -(-x_axis(1,jj)*r(jj,1) - x_axis(2,jj)*r(jj,2) - x_axis(3,jj)*r(jj,3)) /norm(r(jj,:)));
        earth_angle(jj,3) = acosd( -( y_axis(1,jj)*r(jj,1) + y_axis(2,jj)*r(jj,2) + y_axis(3,jj)*r(jj,3)) /norm(r(jj,:)));
        earth_angle(jj,4) = acosd( -(-y_axis(1,jj)*r(jj,1) - y_axis(2,jj)*r(jj,2) - y_axis(3,jj)*r(jj,3)) /norm(r(jj,:)));
        earth_angle(jj,5) = acosd( -( z_axis(1,jj)*r(jj,1) + z_axis(2,jj)*r(jj,2) + z_axis(3,jj)*r(jj,3)) /norm(r(jj,:)));
        earth_angle(jj,6) = acosd( -(-z_axis(1,jj)*r(jj,1) - z_axis(2,jj)*r(jj,2) - z_axis(3,jj)*r(jj,3)) /norm(r(jj,:)));
    end

    figure(16)
    plot(t_vect(:,1)/3600,sun_angle(:,1,1),t_vect(:,1)/3600,sun_angle(:,2,1),t_vect(:,1)/3600,sun_angle(:,3,1),t_vect(:,1)/3600,sun_angle(:,4,1),t_vect(:,1)/3600,sun_angle(:,5,1),t_vect(:,1)/3600,sun_angle(:,6,1))
    hold on
    grid on
    xlabel('Time [h]')
    ylabel('Sun Angles [deg]')
    legend('+x','-x','+y','-y','+z','-z')
    
    figure(17)
    plot(t_vect(:,1)/3600,earth_angle(:,1,1),t_vect(:,1)/3600,earth_angle(:,2,1),t_vect(:,1)/3600,earth_angle(:,3,1),t_vect(:,1)/3600,earth_angle(:,4,1),t_vect(:,1)/3600,earth_angle(:,5,1),t_vect(:,1)/3600,earth_angle(:,6,1))
    hold on
    grid on
    xlabel('Time [h]')
    ylabel('Earth Angles [deg]')
    legend('+x','-x','+y','-y','+z','-z')
end

%% Power Plots
%Panels Power
figure(12)
cmap = jet(1+size(PanelsPower,2));
cmap(end-3,:)=[];
for k = 1:size(PanelsPower,2)
    plot(t_adcs/3600,PanelsPower(1:end,k), 'Color', cmap(k, :));
    hold on
end
xlabel('t [h]')
ylabel('Faces Power [W]')
legend('P_{+x}','P_{-x}','P_{+y}','P_{-y}','P_{+z}','P_{-z}','P_{+w1}','P_{-w1}','P_{+w2}','P_{-w2}')
grid on

%Panels Total Power
figure(13)
plot(t_adcs/3600,PanelsTotalPower(1:end),'k-.',t_adcs/3600,mean(PanelsTotalPower(1:end))*ones(size(t_adcs)),'r')
hold on
plot(t_adcs/3600,mean(PanelsTotalPower(PanelsTotalPower~=0))*ones(size(PanelsTotalPower)),'Color',[0.5,0.8,0])
xlabel('t [h]')
ylabel('Generated Power [W]')
legend('P_{tot}','P_{avg}','P_{avg_{Sun}}')
grid on

%Control Power (Mag)
figure(14)
plot(t_adcs/3600,PMag(1:end,1),t_adcs/3600,PMag(1:end,2),t_adcs/3600,PMag(1:end,3),t_adcs/3600,PMagTotal,'k-.')
xlabel('t [h]')
ylabel('Magnetorquer Power Consumption [W]')
legend('P_{m_1}', 'P_{m_2}', 'P_{m_3}','P_{m_{tot}}')
grid on
hold on

%Control Power (Wheel)
figure(15)
plot(t_adcs/3600,PWheels(1:end,1),t_adcs/3600,PWheels(1:end,2),t_adcs/3600,PWheels(1:end,3),t_adcs/3600,PWheels(1:end,4),t_adcs/3600,PWheelsTotal,'k-.')
xlabel('t [h]')
ylabel('Wheel Power Consumption [W]')
legend('P_{w_1}', 'P_{w_2}', 'P_{w_3}', 'P_{w_4}','P_{w_{tot}}')
grid on
hold on

%% FigTemplate
FigTemplate