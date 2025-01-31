% ----------------------------------------------------------------------------
%
%                           function eci_ecef
%
%  this function transforms between the earth fixed (itrf) frame, and
%    the eci mean equator mean equinox (gcrf).
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    reci        - position vector eci                           km
%    veci        - velocity vector eci                           km/s
%    aeci        - acceleration vector eci                       km/s2
%    direct      - direction                                     eto, efrom
%    iau80arr    - iau76/fk5 eop constants
%    jdtt        - julian date of tt                             days from 4713 bc
%    jdftt       - fractional julian centuries of tt             days
%    jdut1       - julian date of ut1                            days from 4713 bc
%    lod         - excess length of day                          sec
%    xp          - polar motion coefficient                      rad
%    yp          - polar motion coefficient                      rad
%    ddpsi       - delta psi correction to gcrf                  rad
%    ddeps       - delta eps correction to gcrf                  rad
%
%  outputs       :
%    recef       - position vector earth fixed                   km
%    vecef       - velocity vector earth fixed                   km/s
%    aecef       - acceleration vector ecef                      km/s2
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
%    vallado       2022, 223-231
%
% [recef, vecef, aecef] = eci2ecef(reci, veci, aeci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
% ----------------------------------------------------------------------------

function [recef, vecef, aecef] = eci2ecef(reci, veci, aeci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps )
    constastro;

    [fArgs] = fundarg(ttt, '80');

    [prec,psia,wa,ea,xa] = precess ( ttt, '80' );

    [deltapsi, trueeps, meaneps, nut] = nutation  (ttt, ddpsi, ddeps, iau80arr, fArgs);

    [st, stdot] = sidereal(jdut1, deltapsi, meaneps, fArgs(5), lod, eqeterms, '80' );

    [pm] = polarm(xp,yp,ttt,'80');

    thetasa= earthrot * (1.0  - lod/86400.0 );
    omegaearth = [0; 0; thetasa;];

    rpef  = st'*nut'*prec'*reci;
    recef = pm'*rpef;

    vpef  = st'*nut'*prec'*veci - cross( omegaearth,rpef );
    vecef = pm'*vpef;

    temp  = cross(omegaearth,rpef);

    % two additional terms not needed if satellite is not on surface of the Earth
    aecef = pm'*(st'*nut'*prec'*aeci) - cross(omegaearth,temp) - 2.0*cross(omegaearth,vpef);

end