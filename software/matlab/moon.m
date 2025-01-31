% -----------------------------------------------------------------------------
%
%                           function moon
%
%  this function calculates the geocentric equatorial (ijk) position vector
%    for the moon given the julian date.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    jd          - julian date                              days from 4713 bc
%
%  outputs       :
%    rmoon       - ijk position vector of moon              km
%    rtasc       - right ascension                          rad
%    decl        - declination                              rad
%
%  locals        :
%    eclplong    - ecliptic longitude
%    eclplat     - eclpitic latitude
%    hzparal     - horizontal parallax
%    l           - geocentric direction Math.Cosines
%    m           -             "     "
%    n           -             "     "
%    ttdb        - julian centuries of tdb from
%                  jan 1, 2000 12h
%    hr          - hours                                    0 .. 24
%    min         - minutes                                  0 .. 59
%    sec         - seconds                                  0.0  .. 59.99
%    deg         - degrees
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2022, 294, alg 31, ex 5-3
%
% [rmoon, rtasc, decl] = moon ( jd );
% ------------------------------------------------------------------------------

function [rmoon, rtasc, decl] = moon ( jd )

    twopi      =     2.0*pi;
    deg2rad    =     pi/180.0;
    constastro;

    % -------------------------  implementation   -----------------
    ttdb = ( jd - 2451545.0  ) / 36525.0;

    eclplong= 218.32  + 481267.8813 *ttdb ...
        + 6.29 *sin( (134.9 +477198.85 *ttdb)*deg2rad ) ...
        - 1.27 *sin( (259.2 -413335.38 *ttdb)*deg2rad ) ...
        + 0.66 *sin( (235.7 +890534.23 *ttdb)*deg2rad ) ...
        + 0.21 *sin( (269.9 +954397.70 *ttdb)*deg2rad ) ...
        - 0.19 *sin( (357.5 + 35999.05 *ttdb)*deg2rad ) ...
        - 0.11 *sin( (186.6 +966404.05 *ttdb)*deg2rad );      % deg

    eclplat =   5.13 *sin( ( 93.3 +483202.03 *ttdb)*deg2rad ) ...
        + 0.28 *sin( (228.2 +960400.87 *ttdb)*deg2rad ) ...
        - 0.28 *sin( (318.3 +  6003.18 *ttdb)*deg2rad ) ...
        - 0.17 *sin( (217.6 -407332.20 *ttdb)*deg2rad );      % deg

    hzparal =  0.9508  + 0.0518 *cos( (134.9 +477198.85 *ttdb) ...
        *deg2rad ) ...
        + 0.0095 *cos( (259.2 -413335.38 *ttdb)*deg2rad ) ...
        + 0.0078 *cos( (235.7 +890534.23 *ttdb)*deg2rad ) ...
        + 0.0028 *cos( (269.9 +954397.70 *ttdb)*deg2rad );    % deg

    eclplong = rem( eclplong*deg2rad, twopi );
    eclplat  = rem( eclplat*deg2rad, twopi );
    hzparal  = rem( hzparal*deg2rad, twopi );

    obliquity= 23.439291  - 0.0130042 *ttdb;  %deg
    obliquity= obliquity *deg2rad;

    % ------------ find the geocentric direction cosines ----------
    l= cos( eclplat ) * cos( eclplong );
    m= cos(obliquity)*cos(eclplat)*sin(eclplong) ...
        - sin(obliquity)*sin(eclplat);
    n= sin(obliquity)*cos(eclplat)*sin(eclplong) ...
        + cos(obliquity)*sin(eclplat);

    % ------------- calculate moon position vector ----------------
    magr = 1.0 /sin( hzparal );

    rmoon(1)= magr*l;
    rmoon(2)= magr*m;
    rmoon(3)= magr*n;

    % -------------- find rt ascension and declination ------------
    rtasc= atan2( m,l );
    decl = asin( n );

end