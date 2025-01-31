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
% [reci, veci] = tradec2rv(trr, trtasc, tdecl, tdrr, tdrtasc, tddecl, rseci, vseci);
% ------------------------------------------------------------------------------

function [reci, veci] = tradec2rv(trr, trtasc, tdecl, tdrr, tdrtasc, tddecl, rseci, vseci)
    constmath;

    % --------  calculate topocentric slant range vectors ------------------
    rhov(1) = trr * cos(tdecl) * cos(trtasc);
    rhov(2) = trr * cos(tdecl) * sin(trtasc);
    rhov(3) = trr * sin(tdecl);

    drhov(1) = tdrr * cos(tdecl) * cos(trtasc) - ...
        trr * sin(tdecl) * cos(trtasc) * tddecl - ...
        trr * cos(tdecl) * sin(trtasc) * tdrtasc;
    drhov(2) = tdrr * cos(tdecl) * sin(trtasc) - ...
        trr * sin(tdecl) * sin(trtasc) * tddecl + ...
        trr * cos(tdecl) * cos(trtasc) * tdrtasc;
    drhov(3) = tdrr * sin(tdecl) + trr * cos(tdecl) * tddecl;

    % ------ find eci range vector from site to satellite ------
    reci = rhov + rseci;
    veci = drhov + vseci;

end  % tradec2rv
