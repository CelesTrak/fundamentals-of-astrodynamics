%  ------------------------------------------------------------------------------
%
%                           procedure rv_tradec
%
%  this procedure converts topocentric right-ascension declination with
%    position and velocity vectors. the velocity vector is used to find the
%    solution of singular cases.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    reci        - eci position vector                      km
%    veci        - eci velocity vector                      km/s
%    rseci       - eci site position vector                 km
%    direct      - direction to convert                     efrom  eto
%
%  outputs       :
%    rho         - topo radius of the sat                   km
%    trtasc      - topo right ascension                     rad
%    tdecl       - topo declination                         rad
%    drho        - topo radius of the sat rate              km/s
%    tdrtasc     - topo right ascension rate                rad/s
%    tddecl      - topo declination rate                    rad/s
%
%  locals        :
%    rhov        - eci range vector from site               km
%    drhov       - eci velocity vector from site            km/s
%    latgc       - geocentric lat of satellite, not nadir point  -pi/2 to pi/2 rad
%
%  coupling      :
%    mag         - magnitude of a vector
%    addvec      - add two vectors
%    dot         - dot product of two vectors
%
%  references    :
%    vallado       2022, 257, eq 4-1, 4-2, alg 26
%
%  [rho, trtasc, tdecl, drho, dtrtasc, dtdecl] = rv2tradec ( recef, vecef, rsecef, vsecef );
% ------------------------------------------------------------------------------

function [rho, trtasc, tdecl, drho, dtrtasc, dtdecl] = rv2tradec ( recef, vecef, rsecef, vsecef )

    constmath;

    % --------------------- implementation ------------------------
    % ------- find ecef slant range vector from site to satellite ---------
    rhovecef  = recef - rsecef;
    drhovecef = vecef - vsecef;
    rho      = mag(rhovecef);

    % --------------- calculate topocentric rtasc and decl ---------------
    temp = sqrt( rhovecef(1) * rhovecef(1) + rhovecef(2) * rhovecef(2) );
    if (temp < small)
        trtasc = atan2( drhovecef(2), drhovecef(1) );
    else
        trtasc = atan2( rhovecef(2), rhovecef(1) );
    end

    % directly over the north pole
    if (temp < small)
        tdecl = sign(rhovecef(3)) * halfpi;   % +- 90 deg
    else
        magrhoecef = mag(rhovecef);
        tdecl = asin( rhovecef(3) / magrhoecef );
    end
    if (trtasc < 0.0)
        trtasc = trtasc + 2.0*pi;
    end

    % ---------- calculate topcentric rtasc and decl rates -------------
    temp1 = -rhovecef(2) * rhovecef(2) - rhovecef(1) * rhovecef(1);
    drho = dot(rhovecef, drhovecef) / rho;
    if ( abs( temp1 ) > small )
        dtrtasc = ( drhovecef(1)*rhovecef(2) - drhovecef(2) * rhovecef(1) ) / temp1;
    else
        dtrtasc = 0.0;
    end

    if ( abs( temp ) > small )
        dtdecl = ( drhovecef(3) - drho * sin( tdecl ) ) / temp;
    else
        dtdecl = 0.0;
    end

end   % rv2tradec
