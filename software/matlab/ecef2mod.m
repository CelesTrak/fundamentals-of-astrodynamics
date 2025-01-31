% ----------------------------------------------------------------------------
%
%                           function ecef_mod
%
%  this function transforms a vector from earth fixed (itrf) frame to
%  mean of date (mod).
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
%    rmod        - position vector mod                           km
%    vmod        - velocity vector mod                           km/s
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
%    vallado       2022 227
%
%  [rmod, vmod, amod] = ecef2mod( recef, vecef, aecef, iau80arr, ttt,
%       jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps ); 
% ----------------------------------------------------------------------------

function [rmod, vmod, amod] = ecef2mod( recef, vecef, aecef, iau80arr, ...
        ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps )
    constastro;

    [fArgs] = fundarg(ttt, '80');

    % ---- find matrices
    [deltapsi, trueeps, meaneps, nut] = nutation  (ttt, ddpsi, ddeps, iau80arr, fArgs);

    [st, stdot] = sidereal(jdut1, deltapsi, meaneps, fArgs(5), lod, eqeterms, '80' );

    [pm] = polarm(xp,yp,ttt,'80');

    % ---- perform transformations
    thetasa= earthrot * (1.0  - lod/86400.0 );
    omegaearth = [0; 0; thetasa;];

    %trueeps-meaneps
    %deltapsi
    %nut
    rpef = pm*recef;
    rmod = nut*st*rpef;

    vpef = pm*vecef;
    vmod = nut*st*(vpef + cross(omegaearth,rpef));

    temp = cross(omegaearth,rpef);
    amod = nut*st*( pm*aecef + cross(omegaearth,temp) ...
        + 2.0*cross(omegaearth,vpef) );

end