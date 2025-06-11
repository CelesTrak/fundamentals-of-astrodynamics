% ------------------------------------------------------------------------------
%
%                           function sun
%
%  this function calculates the geocentric equatorial position vector
%    the sun given the julian date. Sergey K (2022) has noted that improved results
%    are found assuming the oputput is in a precessing frame (TEME) and converting to ICRF.
%    this is the low precision formula and is valid for years from 1950 to 2050.
%    accuaracy of apparent coordinates is about 0.01 degrees.  notice many of
%    the calculations are performed in degrees, and are not changed until later.
%    this is due to the fact that the almanac uses degrees exclusively in their formulations.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    jd          - julian date  (UTC)                       days from 4713 bc
%
%  outputs       :
%    rsun        - inertial position vector of the sun      au
%    rtasc       - right ascension                          rad
%    decl        - declination                              rad
%
%  locals        :
%    meanlong    - mean longitude
%    meananomaly - mean anomaly
%    eclplong    - ecliptic longitude
%    obliquity   - mean obliquity of the ecliptic
%    tut1        - julian centuries of ut1 from
%                  jan 1, 2000 12h
%    ttdb        - julian centuries of tdb from
%                  jan 1, 2000 12h
%    hr          - hours                                    0 .. 24              10
%    min         - minutes                                  0 .. 59              15
%    sec         - seconds                                  0.0  .. 59.99          30.00
%    temp        - temporary variable
%    deg         - degrees
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2022, 285, alg 29, ex 5-1
%
% [rsun,rtasc,decl] = sun ( jd );
% ------------------------------------------------------------------------------

function [rsun,rtasc,decl] = sun ( jd )

    twopi      =     2.0*pi;
    deg2rad    =     pi/180.0;
    show = 'n';

    % -------------------------  implementation   -----------------
    % -------------------  initialize values   --------------------
    tut1= ( jd - 2451545.0  )/ 36525.0;

    if show == 'y'
        fprintf(1,'tut1 %14.9f \n',tut1);
    end

    meanlong= 280.460  + 36000.771285*tut1;
    meanlong= rem( meanlong,360.0  );  %deg

    ttdb= tut1;
    meananomaly= 357.528  + 35999.0509575 *ttdb;
    meananomaly= rem( meananomaly*deg2rad,twopi );  %rad
    if ( meananomaly < 0.0  )
        meananomaly= twopi + meananomaly;
    end

    eclplong= meanlong + 1.914666471 *sin(meananomaly) ...
        + 0.019994643 *sin(2.0 *meananomaly); %deg
    eclplong= rem( eclplong,360.0  );  %deg

    obliquity= 23.439291  - 0.0130042 *ttdb;  %deg

    eclplong = eclplong *deg2rad;
    obliquity= obliquity *deg2rad;

    % --------- find magnitude of sun vector, )   components ------
    magr= 1.000140612  - 0.016708617 *cos( meananomaly ) ...
        - 0.000139589 *cos( 2.0 *meananomaly );    % in au's

    rsun(1)= magr*cos( eclplong );
    rsun(2)= magr*cos(obliquity)*sin(eclplong);
    rsun(3)= magr*sin(obliquity)*sin(eclplong);

    if show == 'y'
        fprintf(1,'meanlon %11.6f meanan %11.6f eclplon %11.6f obli %11.6f \n', ...
            meanlong,meananomaly/deg2rad,eclplong/deg2rad,obliquity/deg2rad);
        fprintf(1,'rs %11.9f %11.9f %11.9f \n',rsun);
        fprintf(1,'magr %14.7f \n',magr);
    end

    rtasc= atan( cos(obliquity)*tan(eclplong) );

    % --- check that rtasc is in the same quadrant as eclplong ----
    if ( eclplong < 0.0  )
        eclplong= eclplong + twopi;    % make sure it's in 0 to 2pi range
    end
    if ( abs( eclplong-rtasc ) > pi*0.5  )
        rtasc= rtasc + 0.5 *pi*round( (eclplong-rtasc)/(0.5 *pi));
    end
    decl = asin( sin(obliquity)*sin(eclplong) );

end