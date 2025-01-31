% ----------------------------------------------------------------------------
%
%                           function ecef_eci06
%
%  this function transforms between the earth fixed (itrf) frame, and
%    the eci mean equator mean equinox (gcrf).
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    recef       - position vector earth fixed                   km
%    vecef       - velocity vector earth fixed                   km/s
%    aecef       - acceleration vector earth fixed               km/s2
%    enum        - direction                                     eto, efrom
%    iau06arr    - iau2006 eop constants
%    xysarr      - array of xys data records                     rad
%    jdtt        - julian date of tt                             days from 4713 bc
%    ttt         - julian centuries of tt                        centuries
%    jdut1       - julian date of ut1                            days from 4713 bc
%    lod         - excess length of day                          sec
%    xp          - polar motion coefficient                      rad
%    yp          - polar motion coefficient                      rad
%    ddx         - delta x correction to gcrf                    rad
%    ddy         - delta y correction to gcrf                    rad
%
%  outputs       :
%    reci        - position vector eci                           km
%    veci        - velocity vector eci                           km/s
%    aeci        - acceleration vector earth inertial            km/s2
%
%  locals        :
%    eqeterms    - terms for ast calculation                     0,2
%    deltapsi    - nutation angle                                rad
%    trueeps     - true obliquity of the ecliptic                rad
%    meaneps     - mean obliquity of the ecliptic                rad
%    prec        - matrix for mod - eci
%    nut         - matrix for tod - mod
%    st          - matrix for pef - tod
%    stdot       - matrix for pef - tod rate
%    pm          - matrix for ecef - pef
%
%  coupling      :
%   precess      - rotation for precession
%   nutation     - rotation for nutation
%   sidereal     - rotation for sidereal time
%   polarm       - rotation for polar motion
%
%  references    :
%    vallado       2022, 211
%
% [reci, veci, aeci] = ecef2eci06( recef, vecef, aecef, iau06arr, ...
%       xys06table, ttt, jdut1, lod, xp, yp, ddx, ddy, opt1 );
% ----------------------------------------------------------------------------

function [reci, veci, aeci] = ecef2eci06( recef, vecef, aecef, iau06arr, ...
        xys06table, ttt, jdut1, lod, xp, yp, ddx, ddy, opt1 )

    sethelp;
    constastro;
    pnb = eye(3);
    st = eye(3);

    [fArgs06] = fundarg(ttt, '06');

    % ---- ceo based, iau2006
    if ~contains(opt1, 'a') || ~contains(opt1, 'b')
        [x, y, s, pnb] = iau06xys (iau06arr, fArgs06, xys06table, ttt, ddx, ddy, opt1);
        [st, stdot] = sidereal(jdut1, 0.0, 0.0, 0.0, lod, 2, '06' );
    end

    [pm] = polarm(xp, yp, ttt, '06');

    % ---- setup parameters for velocity transformations
    thetasa= earthrot * (1.0  - lod/86400.0 );
    omegaearth = [0; 0; thetasa;];

    % ---- perform transformations
    rpef = pm*recef;
    reci = pnb*st*rpef;

    vpef = pm*vecef;
    veci = pnb*st*(vpef + cross(omegaearth,rpef));

    temp = cross(omegaearth,rpef);
    aeci = pnb*st*( pm*aecef + cross(omegaearth,temp) + 2.0*cross(omegaearth,vpef) );

end