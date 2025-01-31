% ----------------------------------------------------------------------------
%
%
%                           function ecef2llb
%
%  this subroutine converts a geocentric equatorial position vector into
%    latitude and longitude.  geodetic and geocentric latitude are found. the
%    inputs must be ecef. there is no iteration and this is borkowski's method.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r           - ecef position vector                     km
%
%  outputs       :
%    latgc       - geocentric lat of satellite, not nadir point           -pi/2 to pi/2 rad
%    latgd       - geodetic latitude                        -pi/2 to pi/2 rad
%    lon         - longitude (west -)                       0 to 2pi rad
%    hellp       - height above the ellipsoid               km
%
%  locals        :
%    rc          - range of site wrt earth center           er
%    height      - height above earth wrt site              er
%    alpha       - angle from iaxis to point, lst           rad
%    olddelta    - previous value of deltalat               rad
%    deltalat    - diff between delta and
%                  geocentric lat                           rad
%    delta       - declination angle of r in ecef           rad
%    rsqrd       - magnitude of r squared                   er2
%    sintemp     - sine of temp                             rad
%    c           -
%
%  coupling      :
%    mag         - magnitude of a vector
%    gcgd        - converts between geocentric and geodetic latitude
%
%  references    :
%    vallado       2022, 174, alg 12 and alg 13, ex 3-3
%
%
% [latgc,latgd,lon,hellp] = ecef2llb ( r );
% ------------------------------------------------------------------------------

function [latgc,latgd,lon,hellp] = ecef2llb ( r )

    twopi =     2.0*pi;
    small =     0.00000001;

    % -------------------------  implementation   -------------------------
    % ---------------- find longitude value  ----------------------
    temp = sqrt( r(1)*r(1) + r(2)*r(2) );
    if ( abs( temp ) < small )
        rtasc= sign(r(3))*pi*0.5;
    else
        rtasc= atan2( r(2) , r(1) );
    end
    lon  = rtasc;
    if ( abs(lon) >= pi )
        if ( lon < 0.0  )
            lon= twopi + lon;
        else
            lon= lon - twopi;
        end
    end

    a= 6378.1363;
    % b = sign(r(3)) * 6356.75160056;
    b = sign(r(3)) * a * sqrt(1.0 - 0.006694385);  % find semiminor axis of the earth

    % -------------- set up initial latitude value  ---------------
    atemp= 1.0 /(a*temp);
    e= (b*r(3)-a*a+b*b)*atemp;
    f= (b*r(3)+a*a-b*b)*atemp;
    third= 1.0 /3.0;
    p= 4.0 *third*(e*f + 1.0  );
    q= 2.0 *(e*e - f*f);
    d= p*p*p + q*q;

    if ( d > 0.0  )
        nu= (sqrt(d)-q)^third - (sqrt(d)+q)^third;
    else
        sqrtp= sqrt(-p);
        nu= 2.0 *sqrtp*cos( third*acos(q/(p*sqrtp)) );
    end
    g= 0.5 *(sqrt(e*e + nu) + e);
    t= sqrt(g*g + (f-nu*g)/(2.0 *g-e)) - g;

    latgd= atan(a*(1.0 -t*t)/(2.0 *b*t));
    hellp= (temp-a*t)*cos( latgd) + (r(3)-b)*sin(latgd);

    latgc = asin(r(3)/mag(r));   % all locations
    %latgc = gd2gc(latgd);  % surface of the Earth locations

end