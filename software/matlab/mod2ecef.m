%
% ----------------------------------------------------------------------------
%
%                           function mod2ecef
%
%  this function trsnforms a vector from the mean equator of date
%  to an earth fixed (ITRF) frame.  the results take into account
%    the effects of nutation, sidereal time, and polar motion.
%
%  author        : david vallado                  719-573-2600   27 jun 2002
%
%
%  inputs          description                    range / units
%    rmod        - position vector mod            km
%    vmod        - velocity vector mod            km/s
%    amod        - acceleration vector mod        km/s2
%    ttt         - julian centuries of tt         centuries
%    jdut1       - julian date of ut1             days from 4713 bc
%    lod         - excess length of day           sec
%    xp          - polar motion coefficient       arc sec
%    yp          - polar motion coefficient       arc sec
%    eqeterms    - terms for ast calculation      0,2
%    ddpsi       - delta psi correction to gcrf   rad
%    ddeps       - delta eps correction to gcrf   rad
%
%  outputs       :
%    recef       - position vector earth fixed    km
%    vecef       - velocity vector earth fixed    km/s
%    aecef       - acceleration vector earth fixedkm/s2
%
%  locals        :
%    deltapsi    - nutation angle                 rad
%    trueeps     - true obliquity of the ecliptic rad
%    meaneps     - mean obliquity of the ecliptic rad
%    omega       -                                rad
%    nut         - matrix for tod - mod
%    st          - matrix for pef - tod
%    stdot       - matrix for pef - tod rate
%    pm          - matrix for ecef - pef
%
%  coupling      :
%   nutation     - rotation for nutation          mod - tod
%   sidereal     - rotation for sidereal time     tod - pef
%   polarm       - rotation for polar motion      pef - ecef
%
%  references    :
%    vallado       2013, 223-229
%
% [recef, vecef, aecef] = mod2ecef( rmod, vmod, amod, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
% ----------------------------------------------------------------------------

function [recef, vecef, aecef] = mod2ecef( rmod, vmod, amod, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    constastro;

    [fArgs] = fundarg(ttt, '80');

    [deltapsi, trueeps, meaneps, nut] = nutation  (ttt, ddpsi, ddeps, iau80arr, fArgs);

    [st, stdot] = sidereal(jdut1, deltapsi, meaneps, fArgs(5), lod, 2, '80');

    [pm] = polarm(xp,yp,ttt,'80');

    thetasa= earthrot * (1.0  - lod/86400.0 );
    omegaearth = [0; 0; thetasa;];

    rpef  = st'*nut'*rmod;
    recef = pm'*rpef;

    vpef  = st'*nut'*vmod - cross( omegaearth,rpef );
    vecef = pm'*vpef;

    temp  = cross(omegaearth,rpef);

    % two additional terms not needed if satellite is not on surface
    % of the Earth
    aecef = pm'*(st'*nut'*amod)- cross(omegaearth,temp) - 2.0*cross(omegaearth,vpef);

