% ----------------------------------------------------------------------------
%
%                           function ecef_tod
%
%  this function transforms a vector from earth fixed (itrf) frame to
%  true of date (tod).
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    recef       - position vector earth fixed                   km
%    vecef       - velocity vector earth fixed                   km/s
%    direct      - direction of transfer                         eto, efrom
%    iau80arr    - iau80 eop constants
%    ttt         - julian centuries of tt                        centuries
%    jdut1       - julian date of ut1                            days from 4713 bc
%    lod         - excess length of day                          sec
%    xp          - polar motion coefficient                      rad
%    yp          - polar motion coefficient                      rad
%    ddpsi       - delta psi correction to gcrf                  rad
%    ddeps       - delta eps correction to gcrf                  rad
%
%  outputs       :
%    rtod        - position vector eci                           km
%    vtod        - velocity vector eci                           km/s
%
%  locals        :
%    deltapsi    - nutation angle                                rad
%    trueeps     - true obliquity of the ecliptic                rad
%    meaneps     - mean obliquity of the ecliptic                rad
%    nut         - matrix for tod - mod
%    st          - matrix for pef - tod
%    stdot       - matrix for pef - tod rate
%    pm          - matrix for ecef - pef
%
%  coupling      :
%   nutation     - rotation for nutation
%   sidereal     - rotation for sidereal time
%   polarm       - rotation for polar motion
%
%  references    :
%    vallado       2022 226
%
%  [rtod, vtod, atod] = ecef2tod(recef, vecef, aecef, iau80arr, ttt, jdut1, ...
%             lod, xp, yp, eqeterms, ddpsi, ddeps);
% ----------------------------------------------------------------------------

function [rtod, vtod, atod] = ecef2tod(recef, vecef, aecef, iau80arr, ttt, ...
        jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps)
    constastro;

    [fArgs] = fundarg(ttt, '80');

    % ---- find matrices - note nut is only needed for st argument inputs
    [deltapsi, trueeps, meaneps, nut] = nutation  (ttt, ddpsi, ddeps, iau80arr, fArgs);

    [st, stdot] = sidereal(jdut1, deltapsi, meaneps, fArgs(5), lod, eqeterms, '80');

    [pm] = polarm(xp, yp, ttt, '80');

    % ---- perform transformations
    thetasa= earthrot * (1.0  - lod/86400.0 );
    omegaearth = [0; 0; thetasa;];

    rpef = pm*recef;
    rtod = st*rpef;

    vpef = pm*vecef;
    vtod = st*(vpef + cross(omegaearth,rpef));

    temp = cross(omegaearth,rpef);
    atod = st*( pm*aecef + cross(omegaearth,temp) + 2.0*cross(omegaearth,vpef) );

end