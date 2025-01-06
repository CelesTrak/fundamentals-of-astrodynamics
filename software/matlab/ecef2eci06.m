% ----------------------------------------------------------------------------
%
%                           function ecef2eci06
%
%  this function transforms a vector from the earth fixed (itrf) frame, to
%    the eci mean equator mean equinox (gcrf).
%
%  author        : david vallado                  719-573-2600   16 jul 2004
%
%  revisions
%
%  inputs          description                    range / units
%    recef       - position vector earth fixed    km
%    vecef       - velocity vector earth fixed    km/s
%    aecef       - acceleration vector earth fixedkm/s2
%    ttt         - julian centuries of tt         centuries
%    jdut1       - julian date of ut1             days from 4713 bc
%    lod         - excess length of day           sec
%    xp          - polar motion coefficient       arc sec
%    yp          - polar motion coefficient       arc sec
%    option      - which approach to use          a-2000a, b-2000b, c-2000xys
%    ddx         - eop correction for x           rad
%    ddy         - eop correction for y           rad
%
%  outputs       :
%    reci        - position vector eci            km
%    veci        - velocity vector eci            km/s
%    aeci        - acceleration vector eci        km/s2
%
%  locals        :
%    pm          - transformation matrix for itrf-pef
%    st          - transformation matrix for pef-ire
%    nut         - transformation matrix for ire-gcrf
%
%  coupling      :
%   iau00pm      - rotation for polar motion      itrf-pef
%   iau00era     - rotation for earth rotyation   pef-ire
%   iau00xys     - rotation for prec/nut          ire-gcrf

%
%  references    :
%    vallado       2004, 205-219
%
% [reci, veci, aeci] = ecef2eci06( recef, vecef, aecef, iau06arr, xys06table, ttt, jdut1, lod, xp, yp, ddx, ddy, opt1 )
% ----------------------------------------------------------------------------

function [reci, veci, aeci] = ecef2eci06( recef, vecef, aecef, iau06arr, xys06table, ttt, jdut1, lod, xp, yp, ddx, ddy, opt1 )

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



