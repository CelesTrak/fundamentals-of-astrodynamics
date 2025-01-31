% ----------------------------------------------------------------------------
%
%                           function tirs_ecef
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
%    iau06arr    - iau2006 eop constants
%    opt         - method option                                 e80, e96, e06cio
%    ttt         - julian centuries of tt                        centuries
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
%    vallado       2022, 213
%
% [recef, vecef, aecef] = tirs2ecef( rtirs, vtirs, atirs, ttt, xp, yp);
% ----------------------------------------------------------------------------

function [recef, vecef, aecef] = tirs2ecef( rtirs, vtirs, atirs, ttt, xp, yp)
    constastro;

    [pm] = polarm(xp, yp, ttt, '06');

    recef = pm'*rtirs;

    vecef = pm'*vtirs;


    % two additional terms not needed if satellite is not on surface
    % of the Earth
    aecef = pm'*atirs;

end