% ----------------------------------------------------------------------------
%
%                           function eci_ecef06
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
% [recef, vecef, aecef] = eci2ecef06(reci, veci, aeci, iau06arr, xysarr, ttt, jdut1, lod, xp, yp, ddx, ddy, opt1 );
% ----------------------------------------------------------------------------

function [recef, vecef, aecef] = eci2ecef06(reci, veci, aeci, iau06arr, xysarr, ttt, jdut1, lod, xp, yp, ddx, ddy, opt1 )
    constastro;
    %      sethelp;
    pnb = eye(3);
    st = eye(3);

    [fArgs06] = fundarg(ttt, '06');

    % ---- cio based, iau2000
    if not(contains(opt1, 'a')) || not(contains(opt1, 'b'))
        [x, y, s, pnb] = iau06xys (iau06arr, fArgs06, xysarr, ttt, ddx, ddy, opt1);
        [st, stdot] = sidereal(jdut1, 0.0, 0.0, 0.0, 0.0, 0, '06');
    end


    [pm] = polarm(xp, yp, ttt, '06');

    % ---- setup parameters for velocity transformations
    thetasa= earthrot * (1.0  - lod/86400.0 );
    omegaearth = [0; 0; thetasa;];

    rpef  = st'*pnb'*reci;
    recef = pm'*rpef;

    vpef  = st'*pnb'*veci - cross( omegaearth,rpef );
    vecef = pm'*vpef;

    temp  = cross(omegaearth,rpef);
    aecef = pm'*(st'*pnb'*aeci - cross(omegaearth,temp) - 2.0*cross(omegaearth,vpef));

end