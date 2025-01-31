%  ------------------------------------------------------------------------------
%
%                           procedure elatlon_rv
%
%  this procedure converts ecliptic latitude and longitude with position and
%    velocity vectors. uses velocity vector to find the solution of Math.Singular
%    cases.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    rijk        - ijk position vector                      km
%    vijk        - ijk velocity vector                      km / s
%    direction   - which set of vars to output              efrom  eto
%
%  outputs       :
%    rr          - radius of the sat                        km
%    ecllat      - ecliptic latitude                        -PI/2 to PI/2 rad
%    ecllon      - ecliptic longitude                       -2PI to 2PI rad
%    drr         - radius of the sat rate                   km/s
%    decllat     - ecliptic latitude rate                   -PI/2 to PI/2 rad
%    eecllon     - ecliptic longitude rate                  -2PI to 2PI rad
%
%  locals        :
%    obliquity   - obliquity of the ecliptic                rad
%    temp        -
%    temp1       -
%    re          - position vec in eclitpic frame
%    ve          - velocity vec in ecliptic frame
%
%  coupling      :
%    mag         - magnitude of a vector
%    rot1        - rotation about 1st axis
%    dot         - dot product
%    arcsin      - arc Math.Sine function
%    Math.Atan2       - arc tangent function that resolves quadrant ambiguites
%
%  references    :
%    vallado       2022 268, eq 4-15%
% [rijk, vijk] = elatlon2rv ( rr, ecllon, ecllat, drr, decllon, decllat );
%  ------------------------------------------------------------------------------

function [rijk, vijk] = elatlon2rv ( rr, ecllon, ecllat, drr, decllon, decllat )

    % --------------------  implementation   ----------------------
    obliquity= 0.40909280;   %23.439291 /rad

    r(1)= rr*cos(ecllat)*cos(ecllon);
    r(2)= rr*cos(ecllat)*sin(ecllon);
    r(3)= rr*sin(ecllat);

    v(1)= drr*cos(ecllat)*cos(ecllon) ...
        - rr*sin(ecllat)*cos(ecllon)*decllat ...
        - rr*cos(ecllat)*sin(ecllon)*decllon;
    v(2)= drr*cos(ecllat)*sin(ecllon) ...
        - rr*sin(ecllat)*sin(ecllon)*decllat ...
        + rr*cos(ecllat)*cos(ecllon)*decllon;
    v(3)= drr*sin(ecllat) + rr*cos(ecllat)*decllat;

    [rijk] = rot1 ( r, -obliquity );
    [vijk] = rot1 ( v, -obliquity );

end