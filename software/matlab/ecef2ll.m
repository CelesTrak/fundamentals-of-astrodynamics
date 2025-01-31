% ----------------------------------------------------------------------------
%
%
%                           function ecef2ll
%
%  these subroutines convert a geocentric equatorial position vector into
%    latitude and longitude.  geodetic and geocentric latitude are found. the
%    inputs must be ecef.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r           - ecef position vector                     km
%
%  outputs       :
%    latgc       - geocentric lat of satellite, not nadir point  -pi/2 to pi/2 rad
%    latgd       - geodetic latitude                        -pi/2 to pi/2 rad
%    lon         - longitude (west -)                       0 to 2pi rad
%    hellp       - height above the ellipsoid               km
%
%  locals        :
%    temp        - diff between geocentric/
%                  geodetic lat                             rad
%    Math.Sintemp     - Math.Sine of temp                   rad
%    olddelta    - previous value of deltalat               rad
%    rtasc       - right ascension                          rad
%    decl        - declination                              rad
%    i           - index
%
%  coupling      :
%    mag         - magnitude of a vector
%    gcgd        - converts between geocentric and geodetic latitude
%
%  references    :
%    vallado       2022, 174, alg 12 and alg 13, ex 3-3
%
%
% [latgc,latgd,lon,hellp] = ecef2ll ( r );
% ------------------------------------------------------------------------------

function [latgc,latgd,lon,hellp] = ecef2ll ( r )

    twopi      =     2.0*pi;
    small      =     0.00000001;         % small value for tolerances
    re         =     6378.1363;
    eesqrd     =     0.006694385000;     % eccentricity of earth sqrd
    rad = 180.0/pi;

    % -------------------------  implementation   -----------------
    magr = mag( r );

    % ----------------- find longitude value  ---------------------
    temp = sqrt( r(1)*r(1) + r(2)*r(2) );
    if ( abs( temp ) < small )
        rtasc= sign(r(3))*pi*0.5;
    else
        rtasc= atan2( r(2), r(1) );
    end
    lon  = rtasc;
    if ( abs(lon) >= pi )   % mod it ?
        if ( lon < 0.0  )
            lon= twopi + lon;
        else
            lon= lon - twopi;
        end
    end
    decl = asin( r(3) / magr );
    latgd= decl;
    %     fprintf(1,'rd %11.7f rtasc %11.7f lon %11.7f latgd %11.7f \n', temp, rtasc*rad, lon*rad, latgd*rad );

    % ------------- iterate to find geodetic latitude -------------
    i= 1;
    olddelta = latgd + 10.0;

    while ((abs(olddelta-latgd)>=small) && (i<10))
        olddelta= latgd;
        sintemp = sin( latgd );
        c       = re  / (sqrt( 1.0 -eesqrd*sintemp*sintemp ));
        latgd= atan( (r(3)+c*eesqrd*sintemp)/temp );
        % fprintf(1,'%3i  c %11.7f gd %11.7f  \n', i, c, latgd*rad );
        i = i + 1;
    end

    % Calculate height
    if (pi*0.5 - abs(latgd)) < pi/180.0  % 1 deg
        hellp   = (temp/cos(latgd)) - c;
    else
        s = c * (1.0 - eesqrd);
        hellp   = r(3)/sin(latgd) - s;
    end

    latgc = asin(r(3)/magr);   % all locations
    %latgc = gd2gc(latgd);  % surface of the Earth locations

end