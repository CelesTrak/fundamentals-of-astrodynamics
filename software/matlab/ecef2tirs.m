
% ----------------------------------------------------------------------------
%
%                           function ecef2tirs
%
%  this function transforms a vector from the earth fixed itrf frame
%    (itrf), to the tirs frame. the results take into account
%    the effects of  polar motion.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%
%  inputs          description                    range / units
%    recef       - position vector earth fixed    km
%    vecef       - velocity vector earth fixed    km/s
%    aecef       - acceleration vector earth fixedkm/s2
%    xp          - polar motion coefficient       arc sec
%    yp          - polar motion coefficient       arc sec
%    ttt         - julian centuries of tt         centuries
%
%  outputs       :
%    rtirs        - position vector tirs            km
%    vtirs        - velocity vector tirs            km/s
%    atirs        - acceleration vector tirs        km/s2
%
%  locals        :
%
%  coupling      :
%   precess      - rotation for precession        mod - eci
%
%  references    :
%    vallado       2001, 219, eq 3-65 to 3-66
%
% [rtirs, vtirs, atirs] = ecef2tirs  ( recef, vecef, aecef, opt, xp, yp, ttt )
% ----------------------------------------------------------------------------

function [rtirs, vtirs, atirs] = ecef2tirs  ( recef, vecef, aecef, ttt, xp, yp )

    [pm] = polarm(xp, yp, ttt, '06');

    rtirs = pm*recef;

    vtirs = pm*vecef;

    atirs = pm*aecef;


