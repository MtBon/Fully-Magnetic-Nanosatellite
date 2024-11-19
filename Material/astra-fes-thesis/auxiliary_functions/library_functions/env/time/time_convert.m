% ------------------------------------------------------------------------------
%
%                           function convtime
%
%  this function finds the time parameters and julian century values for inputs
%    of utc or ut1. numerous outputs are found as shown in the local variables.
%    because calucations are in utc
%
%  author        : david vallado                  719-573-2600    4 jun 2002
%
%  revisions
%    vallado     - add tcg, tcb, etc                              6 oct 2005
%    vallado     - fix documentation for dut1                     8 oct 2002
%    colagrossi  - simplify for MAT purposes                     22 apr 2020
%
%  inputs          description                    range / units
%    year        - year                           1900 .. 2100
%    mon         - month                          1 .. 12
%    day         - day                            1 .. 28,29,30,31
%    hr          - universal time hour            0 .. 23
%    min         - universal time min             0 .. 59
%    sec         - universal time sec (utc)            0.0  .. 59.999
%    dut1        - delta of ut1 - utc             sec
%    dat         - delta of tai - utc             sec
%
%  outputs       :
%    ut1         - universal time                 sec
%    tut1        - julian centuries of ut1
%    jdut1       - julian date (days only)           days from 4713 bc
%    jdut1Frac   - julian date (fraction of a day)   days from 0 hr of the day
%    utc         - coordinated universal time     sec
%    tai         - atomic time                    sec
%    jdtai       - julian date (days only)           days from 4713 bc
%    jdtaiFrac   - julian date (fraction of a day)   days from 0 hr of the day
%    gps         - gps time                    sec
%    jdgps       - julian date (days only)           days from 4713 bc
%    jdgpsFrac   - julian date (fraction of a day)   days from 0 hr of the day
%    tdt         - terrestrial dynamical time     sec
%    ttdt        - julian centuries of tdt
%    jdtt        - julian date (days only)           days from 4713 bc
%    jdttFrac    - julian date (fraction of a day)   days from 0 hr of the day
%    tdb         - terrestrial barycentric time   sec
%    ttdb        - julian centuries of tdb
%    jdtdb       - julian date of tdb             days from 4713 bc
%    tcb         - celestial barycentric time     sec
%    tcg         - celestial geocentric time      sec
%    jdtdb       - julian date (days only)           days from 4713 bc
%    jdtdbFrac   - julian date (fraction of a day)   days from 0 hr of the day
%
%  locals        :
%    hrtemp      - temporary hours                hr
%    mintemp     - temporary minutes              min
%    sectemp     - temporary seconds              sec
%    localhr     - difference to local time       hr
%    jd          - julian date of request         days from 4713 bc
%    me          - mean anomaly of the earth      rad
%
%
%  references    :
%    vallado       2007, 201, alg 16, ex 3-7
%
% [ut1,tut1,jdut1,jdut1frac,utc,mjdutc,jdutc,jdutcfrac,tai,jdtai,jdtaifrac,gps,jdgps,jdgpsfrac,tt,ttt,jdtt,jdttfrac, ...
%         tdb, ttdb, jdtdb, jdtdbfrac, tcg, jdtcg, jdtcgfrac, tcb, jdtcb, jdtcbfrac ] ...
%   = time_convert ( year, mon, day, hr, min, sec, dut1, dat )
% Validated with vallado exercise 3-7 at p195 - 22/4/2020
% Validated with matlab datetime function and online website http://leapsecond.com/java/gpsclock.htm

% ------------------------------------------------------------------------------

function [ut1,tut1,jdut1,jdut1frac,utc,mjdutc,jdutc,jdutcfrac,tai,jdtai,jdtaifrac,gps,jdgps,jdgpsfrac,tt,ttt,jdtt,jdttfrac, ...
          tdb, ttdb, jdtdb, jdtdbfrac, tcg, jdtcg, jdtcgfrac, tcb, jdtcb, jdtcbfrac ] ...
    = time_convert ( year, mon, day, hr, min, sec, dut1, dat )

% ------------------------  implementation   ------------------
[jdutc, jdutcfrac] =  jday_s(year, mon, day, hr , min, sec);
mjdutc  = jdutc+jdutcfrac - 2400000.5;
%mfme = hr*60.0 + min + sec/60.0;

% ------------------ start if ut1 is known ------------------
utc = hms2sec_s(hr, min, sec );

%%%%%%%%%%%%%%%%%%%%%% ut1
ut1= utc + dut1;
[hrtemp,mintemp,sectemp] = sec2hms_s(ut1);
[jdut1, jdut1frac] = jday_s(year,mon,day, hrtemp, mintemp, sectemp );
tut1= (jdut1+jdut1frac - 2451545.0)/ 36525.0;

%%%%%%%%%%%%%%%%%%%%%% tai
tai= utc + dat;
[hrtemp,mintemp,sectemp] = sec2hms_s(tai);
[jdtai, jdtaifrac] = jday_s( year,mon,day, hrtemp, mintemp, sectemp );

%%%%%%%%%%%%%%%%%%%%%% gps
gps= tai - 19; %sec
[hrtemp,mintemp,sectemp] = sec2hms_s(gps);
[jdgps, jdgpsfrac] = jday_s( year,mon,day, hrtemp, mintemp, sectemp );

%%%%%%%%%%%%%%%%%%%%%% tt
tt= tai + 32.184;   % sec
[hrtemp,mintemp,sectemp] = sec2hms_s( tt );
[jdtt, jdttfrac] = jday_s( year,mon,day, hrtemp, mintemp, sectemp);
ttt= (jdtt+jdttfrac - 2451545.0  )/ 36525.0;

%%%%%%%%%%%%%%%%%%%%%% tdb
% usno circular approach
tdb = tt + 0.001657*sin(628.3076*ttt+6.2401) ...
    + 0.000022*sin(575.3385*ttt+4.2970) ...
    + 0.000014*sin(1256.6152*ttt+6.1969) ...
    + 0.000005*sin(606.9777*ttt+4.0212) ...
    + 0.000005*sin(52.9691*ttt+0.4444) ...
    + 0.000002*sin(21.3299*ttt+5.5431) ...
    + 0.000010*ttt*sin(628.3076*ttt+4.2490);  % USNO circ (14)
[hrtemp,mintemp,sectemp] = sec2hms_s( tdb );
[jdtdb, jdtdbfrac] = jday_s( year,mon,day, hrtemp, mintemp, sectemp );
ttdb = (jdtdb+jdtdbfrac - 2451545.0  )/ 36525.0;

%%%%%%%%%%%%%%%%%%%%%% tcg
% approx with tai
tcg = tt + 6.969290134e-10*(jdtai - 2443144.5)*86400.0;  % AAS 05-352 (10) and IERS TN (104)
[hrtemp,mintemp,sectemp] = sec2hms_s( tcg );
[jdtcg, jdtcgfrac] = jday_s( year,mon,day, hrtemp, mintemp, sectemp );
%tt2 = tcg-6.969290134e-10*(jdtcg+jdtcgfrac-2443144.5003725)*86400.0;

%%%%%%%%%%%%%%%%%%%%%% tcb
tcbmtdb = 1.55051976772e-8*(jdtai+jdtaifrac - 2443144.5)*86400.0;  % sec, value for de405 AAS 05-352 (10) and IERS TN (104)
tcb = tdb + tcbmtdb;
[hrtemp,mintemp,sectemp] = sec2hms_s( tcb );
[jdtcb, jdtcbfrac] = jday_s( year,mon,day, hrtemp, mintemp, sectemp );
end



% -----------------------------------------------------------------------------
%
%                           function jday.m
%
%  this function finds the julian date given the year, month, day, and time.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    year        - year                           1900 .. 2100
%    mon         - month                          1 .. 12
%    day         - day                            1 .. 28,29,30,31
%    hr          - universal time hour            0 .. 23
%    min         - universal time min             0 .. 59
%    sec         - universal time sec             0.0 .. 59.999
%    whichtype   - julian .or. gregorian calender   'j' .or. 'g'
%
%  outputs       :
%    jd          - julian date                    days from 4713 bc
%    jdfrac      - julian date fraction of a day   0.0 to 1.0
%
%  locals        :
%    none.
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2007, 189, alg 14, ex 3-14
%
% [jd, jdfrac] = jday(yr, mon, day, hr, min, sec)
% -----------------------------------------------------------------------------

function [jd, jdfrac] = jday_s(yr, mon, day, hr, min, sec)

% ------------------------  implementation   ------------------
jd = 367.0 * yr  ...
    - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 )   ...
    + floor( 275 * mon / 9.0 ) ...
    + day + 1721013.5;   % use - 678987.0 to go to mjd directly
jdfrac = (sec + min * 60.0 + hr *3600.0) / 86400.0;

% check jdfrac
if jdfrac > 1.0
    jd = jd + floor(jdfrac);
    jdfrac = jdfrac - floor(jdfrac);
end

%  - 0.5 * sign(100.0 * yr + mon - 190002.5) + 0.5;
end

% -----------------------------------------------------------------------------
%
%                           function sec2hms
%
function [hr,min,sec] = sec2hms_s(utsec)

% ------------------------  implementation   ------------------
temp  = utsec / 3600.0;
hr    = fix( temp );
min   = fix( (temp - hr)* 60.0 );
sec   = (temp - hr - min/60.0 ) * 3600.0;
end

% -----------------------------------------------------------------------------
%
%                           function hms2sec
%
function [utsec] = hms2sec_s(hr,min,sec )

% ------------------------  implementation   ------------------
utsec  = hr * 3600.0 + min * 60.0 + sec;

end
