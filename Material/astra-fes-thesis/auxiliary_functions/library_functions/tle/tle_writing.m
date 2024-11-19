%                            procedure TLE_WRITING
%
% This script generates a TLE file for a given set of data and satellite
%
%
% (c) 2020 MAT @ PoliMi
%
%  ----------------------------------------------------------------------------*/
clc
clear

% Initialize Variables
re = 6378.137;         % km
mu = 398600.4418; %km3/s2 Standard gravitational parameter for the earth

%% TLE File
tlefileName='HERMES_FM01.txt';

%% TLE Data
% Title
tleTitle='HERMES (FM01)';
% Satellite Number
satID=10069;
% Classification (U=unclassified)
classification='U';
% Launch Year
launchYr=2020;
% Launch number of the year
launchNo=56;
% Piece of the launch
launchPiece='A';
% Ephemeris type (Zero for distributed TLE)
ephTyp=0;
% Tle number
tleN=1;
% Revolution Number at Epoch
revN=0;

%% TLE Epoch
Year=2020; Month=4; Day=1; Hr=12; Min=30; Sec=1;

%% TLE delta Epoch and Number of TLE
NT=24; %Number of TLE to Produce
DT=1/24; %[d] Fraction of Day

%% Spacecraft Ballistic Coefficent (Inverse): CD*A/m [m2/kg]
A=((0.1*0.3)+((0.4+0.1*sqrt(2))*0.3))/2; %Mean S/C Area
BC=2.2*A/6.42;

%% Nominal Keplerian Parameters
a=6929; %km
e=0.0017; %[-]
i=1e-4; %[deg]
W=0; %[deg]
w=90; %[deg]
ta=270.0627; %[deg]

%% Uncertainties of Keplerian Parameters (3sigm)
aU=0*1; %km
eU=0*1e-7; %[-]
iU=0*1e-4; %[deg]
WU=0*1e-4; %[deg]
wU=0*1e-4; %[deg]
taU=0*1e-4; %[deg]

%  ----------------------------------------------------------------------------*/
%  ---------------------------- TLE LOOP ---------------------------------------
%  ----------------------------------------------------------------------------*/

%% Initialize Loop Data
% Keplerian Parameters
a=a+aU*randn/3; %km
e=max(e+eU*randn/3,0); %[-]
i=max(i+iU*randn/3,0); %[deg]
W=wrapTo360(W+WU*randn/3); %[deg]
w=wrapTo360(w+wU*randn/3); %[deg]
ta=wrapTo360(ta+taU*randn/3); %[deg]

fileID = fopen(tlefileName,'w');
for j=1:NT
    
    %  ----------------------------------------------------------------------------*/
    %  ------------------------ Process Data ---------------------------------------
    %  ----------------------------------------------------------------------------*/
    %% Orbit Data
    % Mean Motion
    n=sqrt(mu/a^3); %rad/s
    
    % Bstar from Vallado: Bstar=1/2*(CD*A/m)*rho_0=1/2*BC*rho_0 (adimensionale)
    % Bstar=1/2*BC*(2.461e-5*6378.135)=BC/12.7416
    % From Vallado - Page 140 (eq 2-61)
    Bstar=BC/12.741621; % [-]
    
    % First Derivative of mean motion
    % Assume Constant Density rho_0 = 5E-13 kg/m3
    ndot=3*n/(2*(a*1000))*(5e-13*BC*((mu*1e9)/(a*1000))*(sqrt(1+e^2+2*e*cosd(ta))/(n*sqrt(1-e)))); %rad/s2
    
    % Second derivative of Mean Motion
    % Assume Delta T = 1s - Assume M=TA (e<<<1)
    nk=n;
    ak=a-(2*a)/(3*n)*ndot*1;
    ek=e-(2*(1-e))/(3*n)*ndot*1;
    ndotk=3*nk/(2*(ak*1000))*(5e-13*BC*((mu*1e9)/(ak*1000))*(sqrt(1+ek^2+2*ek*cosd(ta+rad2deg(nk)*1))/(nk*sqrt(1-ek)))); %rad/s2
    nddot=(ndotk-ndot)/1; %rad/s3
    
    %  ----------------------------------------------------------------------------*/
    %  -------------------- Compose TLE strings ------------------------------------
    %  ----------------------------------------------------------------------------*/
    %% Title Data
    tleTitle = pad(tleTitle,69);
    
    %% First Line Data
    
    % Satellite Number - 5 digits
    satIDS=num2str(satID,'%05i');
    % Launch Year - Last 2 digits
    launchYrS=num2str(launchYr,'%02i');
    launchYrS=launchYrS(end-1:end);
    % Launch number of the year - 3 digits
    launchNoS=num2str(launchNo,'%03i');
    % Piece of the launch - up to 3 letters
    launchPieceS = pad(launchPiece,3);
    % Internation Designator
    intDsgS=[launchYrS,launchNoS,launchPieceS];
    % Ephemeris type - 1 digit (Zero for distributed TLE)
    ephTypS=num2str(ephTyp,'%01i');
    % Tle number - 3 digits
    tleNS=num2str(tleN,'%04i');
    
    
    % Epoch
    EPOCH=datetime(Year,Month,Day,Hr,Min,Sec);
    % Epoch year - Last 2 digits
    epYrS=num2str(year(EPOCH),'%02i');
    epYrS=epYrS(end-1:end);
    % Epoch - Day of the year and fractional portion of the day - 3.8 format (12 digits)
    % Day of year - 3 digits
    epD=day(EPOCH,'dayofyear');
    % Fraction of Day
    epF=(hour(EPOCH)/24)+(minute(EPOCH)/24/60)+(second(EPOCH)/24/60/60);
    % Epoch Day
    epDFS=num2str(epD+epF,'%012.8f');
    
    % First derivative of Mean Motion - sign.8 format (10 digits)
    % ndot/2 in rev/day2
    meanMdot=(ndot*86400*86400/(2*pi))/2;
    meanMdotS=num2str(meanMdot,'%+010.8f');
    meanMdotS(2)='';
    
    % Second derivative of Mean Motion - sign 5 digits + 2 digits for sign and exponent (8 digits)
    % nddot/6 in rev/day3
    meanMddot=(nddot*86400*86400*86400/(2*pi))/6;
    % Keep just one digit exponent in this number
    if abs(meanMddot)<1e-5
        meanMddotS='+00000-0';
    else
        expmeanMddot=floor(log10(meanMddot))+1;
        basemeanMddot=meanMddot/(10^(expmeanMddot));
        % Convert to string
        expmeanMddot=num2str(expmeanMddot,'%+02i');
        basemeanMddot=num2str(basemeanMddot,'%+07.5f');
        basemeanMddot(2:3)='';
        meanMddotS=[basemeanMddot,expmeanMddot];
    end
    
    % Drag term Bstar - sign 5 digits + 2 digits for sign and exponent (8 digits)
    % Keep just one digit exponent in this number
    if abs(Bstar)<1e-9
        BstarS='+00000-0';
    else
        expBstar=floor(log10(Bstar))+1;
        baseBstar=Bstar/(10^(expBstar));
        % Convert to string
        expBstar=num2str(expBstar,'%+02i');
        baseBstar=num2str(baseBstar,'%+07.5f');
        baseBstar(2:3)='';
        BstarS=[baseBstar,expBstar];
    end
    
    %% Second Line Data
    
    % Inclination [deg] 3.4 format (8 digits)
    incS=num2str(i,'%08.4f');
    
    % RAAN [deg] 3.4 format (8 digits)
    raanS=num2str(W,'%08.4f');
    
    % ECC [-] decimal points (7 digits)
    eccS=num2str(e,'%08.7f');
    eccS(1:2)='';
    
    % Argument of Perigee [deg] 3.4 format (8 digits)
    argperS=num2str(w,'%08.4f');
    
    % Mean Anomaly [deg] 3.4 format (8 digits)
    sine= ( sqrt( 1.0 - e*e ) * sind(ta) ) / ( 1.0 + e*cosd(ta) );
    cose= ( e + cosd(ta) ) / ( 1.0  + e*cosd(ta) );
    e0  = atan2d( sine,cose );
    ma  = wrapTo360(e0 - e*sind(e0));
    maS=num2str(ma,'%08.4f');
    
    % Mean Motion [rev/d] 2.8 format (11 digits)
    meanM=n*86400/(2*pi);
    meanMS=num2str(meanM,'%011.8f');
    
    % Revolution Number at Epoch (5 digits)
    revNS=num2str(floor(revN),'%05i');
    
    
    %% Compose String
    TLE0=tleTitle;
    
    TLE1t=['1 ',satIDS,classification,' ',intDsgS,' ',epYrS,epDFS,' ',meanMdotS,' ',meanMddotS,' ',BstarS,' ',ephTypS,' ',tleNS];
    TLE1=[TLE1t,chksum(TLE1t)];
    
    TLE2t=['2 ',satIDS,' ',incS,' ',raanS,' ',eccS,' ',argperS,' ',maS,' ',meanMS,revNS];
    TLE2=[TLE2t,chksum(TLE2t)];
    
    %% Save to file
    fprintf(fileID,'%s\n',TLE0);
    fprintf(fileID,'%s\n',TLE1);
    fprintf(fileID,'%s\n',TLE2);
    
    %% Update Cycle
    % Epoch
    EPOCH=EPOCH+days(DT);
    Year=year(EPOCH);
    Month=month(EPOCH);
    Day=day(EPOCH);
    Hr=hour(EPOCH);
    Min=minute(EPOCH);
    Sec=second(EPOCH);
    % Number of TLE
    tleN=tleN+1;
    % Number of revolutions
    revN=revN+meanM*DT;
    % Propagate Orbit
    [a,e,i,W,w,ta] = TLEpropagation(a,e,i,W,w,ma,DT*86400,ndot,nddot);
    
end
fclose(fileID);



%% TLE Check Function
% Checksum (Modulo 10)
% Letters, blanks, periods, plus signs = 0; minus signs = 1
function cs = chksum(str)
c = 0;
for k = 1:68
    if str(k) > '0' && str(k) <= '9'
        c = c + str(k) - 48;
    elseif str(k) == '-'
        c = c + 1;
    end
end
cs=char(mod(c,10)+48);
end

%% Propagation Function
function [a,ecc,incl,omega,argp,nu] = TLEpropagation(a,ecc,incl,omega,argp,m, dtsec, ndot,nddot)
% KP [km] and [deg]
% dtsec [s]
% ndot and nddot [rad/s^x]

% Constants
re = 6378.137;   %km
mu = 398600.4418; %km3/s2 Standard gravitational parameter for the earth
j2 =  0.00108263; % -

% Compute Derived Quantities
% ndot and nddot are the real value in rad and s - NOT those in TLE divided ny 2 or 6
n= sqrt(mu/(a*a*a)); %rad/s
p= a*(1.0 - ecc*ecc); %km

% ------------- find the value of j2 perturbations -------------
j2op2   = (n*1.5*re^2*j2) / (p*p); %rad/s
omegadot= -j2op2 * cosd(incl); %rad/s
argpdot =  j2op2 * (2.0-2.5*sind(incl)*sind(incl)); %rad/s
mdot    =  n; %rad/s

a     = a - 2.0*ndot*dtsec * a / (3.0*n); %km
ecc   = ecc - 2.0*(1.0 - ecc)*ndot*dtsec / (3.0*n); % -
ecc   = max(ecc,0); % -
%p     = a*(1.0 - ecc*ecc); %km

% ----- update the orbital elements for each orbit type --------
% --- elliptical, parabolic, hyperbolic inclined --
omega= omega + rad2deg(omegadot) * dtsec; %deg
omega= wrapTo360(omega); %deg
argp = argp  + rad2deg(argpdot)  * dtsec; %deg
argp = wrapTo360(argp); %deg
m    = m + rad2deg(mdot)*dtsec + rad2deg(ndot/2)*dtsec*dtsec + rad2deg(nddot/6)*dtsec*dtsec*dtsec; %deg
m    = wrapTo360(m); %deg
% -------------------  Solve Kepler's   -----------------
% Algorithm in [rad]
% Convert mean anomaly in [rad]
m=deg2rad(m); %rad
% -------------------- Elliptical ----------------------
if ( ecc > 1e-8 )
    if ( ((m < 0.0 ) && (m > -pi)) || (m > pi) )
        e0= m - ecc; %rad
    else
        e0= m + ecc; %rad
    end
    ktr= 1;
    e1 = e0 + ( m - e0 + ecc*sin(e0) ) / ( 1.0  - ecc*cos(e0) ); %rad
    while (( abs(e1-e0) > 1e-8 ) && ( ktr <= 50 ))
        ktr = ktr + 1;
        e0= e1; %rad
        e1= e0 + ( m - e0 + ecc*sin(e0) ) / ( 1.0  - ecc*cos(e0) ); %rad
    end
    % -------------  find true anomaly  ---------------
    sinv= ( sqrt( 1.0 -ecc*ecc ) * sin(e1) ) / ( 1.0 -ecc*cos(e1) );
    cosv= ( cos(e1)-ecc ) / ( 1.0  - ecc*cos(e1) );
    nu  = atan2( sinv,cosv ); %rad
else
    % -------------------- circular -------------------
    nu= m; %rad
end
% Convert True and Mean Anomaly in [deg]
nu=rad2deg(nu); %deg

end














