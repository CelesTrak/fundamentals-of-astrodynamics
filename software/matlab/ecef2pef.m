% ----------------------------------------------------------------------------
%
%                           function ecef_pef
%
%  this function transforms a vector from earth fixed (itrf) frame to
%  pseudo earth fixed (pef).
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    recef       - position vector earth fixed                   km
%    vecef       - velocity vector earth fixed                   km/s
%    direct      - direction of transfer                         eto, efrom
%    opt         - method option                                 e80, e96, e06cio
%    ttt         - julian centuries of tt                        centuries
%    lod         - excess length of day                          sec
%    xp          - polar motion coefficient                      rad
%    yp          - polar motion coefficient                      rad
%
%  outputs       :
%    rpef        - position vector pef                           km
%    vpef        - velocity vector pef                           km/s
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
%   polarm       - rotation for polar motion
%
%  references    :
%    vallado       2022, 224
%
% [rpef,vpef,apef] = ecef2pef  ( recef,vecef,aecef, opt, xp, yp, ttt );
% ----------------------------------------------------------------------------

function [rpef,vpef,apef] = ecef2pef  ( recef, vecef, aecef, opt, ttt, xp, yp )

    [pm] = polarm(xp, yp, ttt, opt);

    rpef = pm*recef;

    vpef = pm*vecef;

    apef = pm*aecef;

end