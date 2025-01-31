% ----------------------------------------------------------------------------
%
%                           function eci_pef
%
%  this function transforms between the eci mean equator mean equinox (gcrf), and
%    the pseudo earth fixed frame (pef).
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    reci        - position vector eci                           km
%    veci        - velocity vector eci                           km/s
%    direct      - direction                                     eto, efrom
%    iau80arr    - iau76/fk5 eop constants
%    jdtt        - julian date of tt                             days from 4713 bc
%    jdftt       - fractional julian centuries of tt             days
%    jdut1       - julian date of ut1                            days from 4713 bc
%    lod         - excess length of day                          sec
%    ddpsi       - delta psi correction to gcrf                  rad
%    ddeps       - delta eps correction to gcrf                  rad
%
%  outputs       :
%    rpef       - position vector pef                            km
%    vpef       - velocity vector pef                            km/s
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
%
%  coupling      :
%   precess      - rotation for precession
%   nutation     - rotation for nutation
%   sidereal     - rotation for sidereal time
%
%  references    :
%    vallado       2022, 224
%
% [reci, veci, aeci] = pef2eci(rpef, vpef, apef, iau80arr, ttt, jdut1, lod, eqeterms, ddpsi, ddeps)
% ----------------------------------------------------------------------------

function [reci, veci, aeci] = pef2eci(rpef, vpef, apef, iau80arr, ttt, jdut1, lod, eqeterms, ddpsi, ddeps)
    constastro;
    pnb = zeros(3,3);
    st = zeros(3,3);

    [fArgs] = fundarg(ttt, '80');

    [prec, psia, wa, ea, xa] = precess ( ttt, '80' );

    [deltapsi, trueeps, meaneps, nut] = nutation  (ttt, ddpsi, ddeps, iau80arr, fArgs);
    [st, stdot] = sidereal(jdut1, deltapsi, meaneps, fArgs(5), lod, eqeterms, '80' );

    thetasa= earthrot * (1.0  - lod/86400.0 );
    omegaearth = [0; 0; thetasa;];

    reci = prec*nut*st*rpef;

    veci = prec*nut*st*(vpef + cross(omegaearth,rpef));

    temp = cross(omegaearth,rpef);
    aeci = prec*nut*st*(apef + cross(omegaearth,temp) ...
        + 2.0*cross(omegaearth,vpef));

end