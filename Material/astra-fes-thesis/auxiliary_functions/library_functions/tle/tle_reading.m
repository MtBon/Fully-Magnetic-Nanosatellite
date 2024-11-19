%                            procedure TLE_Reading
%
% This function converts the two line element set character string data to
%    variables. several intermediate varaibles
%    and quantities are determined.
%
% TLE structure contains all the TLE processed data
%
% Error codes
%  0 = No error
% -1 = First line does not correspond
% -2 = Second line does not correspond
% -3 = Satnum does not match
% -4 = Checksum error in first line
% -5 = Checksum error in second line
%
% (c) 2020 MAT @ PoliMi
%
%  ----------------------------------------------------------------------------*/

% Initialize Variables
mu = 3.986004362330000e+14; %  Standard gravitational parameter for the earth %m3/s2
deg2rad  =   pi / 180.0;         %  0.01745329251994330;  % [deg/rad]

% File Name
tlefileName='HERMES_FM01';
tleID = fopen([tlefileName,'.txt'],'r');

% Initialize Loop
id=1;
processing=1;
TLEstr0 = string(fgetl(tleID));
TLEstr1 = string(fgetl(tleID));
TLEstr2 = string(fgetl(tleID));
while processing
    error=zeros(1,5);
    
    % Convert to Char Arrays
    TLEstr1=char(TLEstr1);
    TLEstr2=char(TLEstr2);
    
    % Check first line
    if ~chksum(TLEstr1)
        error(4)=-4;
        warning('Checksum of first line fails at TLE n. %i',id);
    end
    % Check second line
    if ~chksum(TLEstr2)
        error(5)=-5;
        warning('Checksum of second line fails at TLE n. %i',id);
    end
    
    %// set the implied decimal points since doing a formated read
    %// fixes for bad input data values (missing, ...)
    for j = 11:15
        if (TLEstr1(j) == ' ')
            TLEstr1(j) = '_';
        end
    end
    
    if (TLEstr1(45) ~= ' ')
        TLEstr1(44) = TLEstr1(45);
    end
    TLEstr1(45) = '.';
    
    if (TLEstr1(8) == ' ')
        TLEstr1(8) = 'U';
    end
    
    if (TLEstr1(10) == ' ')
        TLEstr1(10) = '.';
    end
    
    for j = 46:50
        if (TLEstr1(j) == ' ')
            TLEstr1(j) = '0';
        end
    end
    
    if (TLEstr1(52) == ' ')
        TLEstr1(52) = '0';
    end
    
    if (TLEstr1(54) ~= ' ')
        TLEstr1(53) = TLEstr1(54);
    end
    
    TLEstr1(54) = '.';
    
    TLEstr2(26) = '.';
    
    for j = 27:33
        if (TLEstr2(j) == ' ')
            TLEstr2(j) = '0';
        end
    end
    
    if (TLEstr1(63) == ' ')
        TLEstr1(63) = '0';
    end
    
    if ((length(TLEstr1) < 68) || (TLEstr1(68) == ' '))
        TLEstr1(68) = '0';
    end
    
    % Initialize tle variables
    tle.date=(zeros(1,6));
    tle.julian=(zeros(1,1));
    tle.elements=(zeros(1,9));
    tle.error=(zeros(1,5));
    
    % parse first line
    cardnumb1 = real(str2double(TLEstr1(1)));
    if cardnumb1~=1
        error(1)=-1;
        warning('First line ID error at TLE n. %i',id);
    end
    tle.satnum1 = real(str2double(TLEstr1(3:7)));
    tle.classification = (TLEstr1(8));
    tle.intldesg = (TLEstr1(10:17));
    tle.epochyr = real(str2double(TLEstr1(19:20)));
    tle.epochdays = real(str2double(TLEstr1(21:32)));
    tle.ndot = real(str2double(TLEstr1(34:43)));
    tle.nddot = real(str2double(TLEstr1(44:50)));
    tle.nexp = real(str2double(TLEstr1(51:52)));
    tle.bstar = real(str2double(TLEstr1(53:59)));
    tle.ibexp = real(str2double(TLEstr1(60:61)));
    tle.numb = real(str2double(TLEstr1(63)));
    tle.elnum = real(str2double(TLEstr1(65:68)));
    tle.checksum1 = real(str2double(TLEstr1(69)));
    
    % parse second line
    cardnumb2 = real(str2double(TLEstr2(1)));
    if cardnumb2~=2
        error(2)=-2;
        warning('Second line ID error at TLE n. %i',id);
    end
    tle.satnum2 = real(str2double(TLEstr2(3:7)));
    tle.inclo = real(str2double(TLEstr2(8:16)));
    tle.nodeo = real(str2double(TLEstr2(17:25)));
    tle.ecco = real(str2double(TLEstr2(26:33)));
    tle.argpo = real(str2double(TLEstr2(34:42)));
    tle.mo = real(str2double(TLEstr2(43:51)));
    tle.no_rpd = real(str2double(TLEstr2(52:63)));
    tle.revnum = real(str2double(TLEstr2(64:68)));
    tle.checksum2 = real(str2double(TLEstr2(69)));
    
    % check satnum
    if tle.satnum1~=tle.satnum2
        error(3)=-3;
        warning('Satellite ID error at TLE n. %i',id);
    end
    
    %// ---- find no, ndot, nddot ----
    tle.nddot= tle.nddot * 10.0^tle.nexp;
    tle.bstar= tle.bstar * 10.0^tle.ibexp;
    
    n     = (tle.no_rpd * (2.0*pi) / (24.0*60.0*60.0)); % Mean motion %rad/s
    ndot  = (tle.ndot * (2.0*pi) / (24.0*60.0*60.0)^2); % ndot/2 - 1st Derivative Mean motion %rad/s2
    nddot = (tle.nddot * (2.0*pi) / (24.0*60.0*60.0)^3); % nddot/6 - 2nd Derivative Mean motion %rad/s3
    
    %// ------- keplerian parameter ----
    
    a  = ((mu/n^2)^(1/3));     % Semi-major axis [m]
    e  = (tle.ecco); % Eccentricity [-]
    i  = (tle.inclo  * deg2rad); % Inclination %rad
    w  = (tle.argpo  * deg2rad); % Argument of pericenter %rad
    OM = (tle.nodeo * deg2rad); % Right Ascension of Ascending node %rad
    M  = (tle.mo * deg2rad); % Mean anomaly %rad
    
    %// ------------- temp fix for years from 1957-2056 ----------------
    %// ------ correct fix will occur when year is 4-digit in 2le ------
    if (tle.epochyr < 57)
        year= tle.epochyr + 2000;
    else
        year= tle.epochyr + 1900;
    end
    
    % --------------- set up array of days in month  --------------
    lmonth=([31,28,31,30,31,30,31,31,30,31,30,31]);
    if rem(year-1900,4) == 0
        lmonth(2)= (29);
    end
    
    % ----------------- find month and day of month ---------------
    dayofyr= floor(tle.epochdays);
    j= 1;
    inttemp= 0;
    while ( dayofyr > inttemp + lmonth(j) ) && ( j < 12 )
        inttemp= inttemp + lmonth(j);
        j = j+1;
    end
    mon= j;
    day= dayofyr - inttemp;
    
    % ----------------- find hours minutes and seconds ------------
    temp = (tle.epochdays - dayofyr )*24.0;
    hr  = fix( temp );
    temp = (temp-hr) * 60.0;
    minute = fix( temp );
    sec = (temp-minute) * 60.0;
    
    % ------------------  compute julian day epoch  ------------------
    jd = 367.0 * year  ...
        - floor( (7 * (year + floor( (mon + 9) / 12.0) ) ) * 0.25 )   ...
        + floor( 275 * mon / 9.0 ) ...
        + day + 1721013.5;   % use - 678987.0 to go to mjd directly
    jdfrac = (sec + minute * 60.0 + hr *3600.0) / 86400.0;
    % check jdfrac
    if jdfrac > 1.0
        jd = jd + floor(jdfrac);
        jdfrac = jdfrac - floor(jdfrac);
    end
    tle.date=([year,mon,day,hr,minute,sec]);
    tle.julian = jd + jdfrac;
    tle.elements=[a,e,i,w,OM,M,n,ndot,nddot];
    tle.error=error;
    % Stop processing if error is found
    if any(error)
        processing=0;
    end
    
    % SAVE TLE
    TLE{id}=tle; %#ok<SAGROW>
    clear tle
    
    % Update Loop
    id=id+1;
    TLEstr0 = string(fgetl(tleID));
    if strcmp(TLEstr0,"-1")
        processing=0;
    else
        TLEstr1 = string(fgetl(tleID));
        TLEstr2 = string(fgetl(tleID));
    end
end
fclose(tleID);
%Save to MAT
save(tlefileName,'TLE')


%% TLE Check Function
% Checksum (Modulo 10)
% Letters, blanks, periods, plus signs = 0; minus signs = 1
function result = chksum(str)
result = false; c = 0;

for k = 1:68
    if str(k) > '0' && str(k) <= '9'
        c = c + str(k) - 48;
    elseif str(k) == '-'
        c = c + 1;
    end
end

if mod(c,10) == str(69) - 48
    result = true;
end

end





























