% ----------------------------------------------------------------------------
%
%                           function eci2ecef06
%
%  this function trsnforms a vector from the mean equator mean equniox frame
%    (gcrf), to an earth fixed (itrf) frame.  the results take into account
%    the effects of precession, nutation, sidereal time, and polar motion.
%
%  author        : david vallado                  719-573-2600   16 jul 2004
%
%  revisions
%
%  inputs          description                    range / units
%    reci        - position vector eci            km
%    veci        - velocity vector eci            km/s
%    aeci        - acceleration vector eci        km/s2
%    ttt         - julian centuries of tt         centuries
%    jdut1       - julian date of ut1             days from 4713 bc
%    lod         - excess length of day           sec
%    xp          - polar motion coefficient       rad
%    yp          - polar motion coefficient       rad
%    option      - which approach to use          a-2000a, b-2000b, c-2000xys
%    ddx         - eop correction for x           rad
%    ddy         - eop correction for y           rad
%
%  outputs       :
%    recef       - position vector earth fixed    km
%    vecef       - velocity vector earth fixed    km/s
%    aecef       - acceleration vector earth fixedkm/s2
%
%  locals        :
%    pm          - transformation matrix for itrf-pef
%    st          - transformation matrix for pef-ire
%    nut         - transformation matrix for ire-gcrf
%
%  coupling      :
%   iau00pm      - rotation for polar motion      itrf-pef
%   iau00era     - rotation for earth rotation    pef-ire
%   iau00xys     - rotation for prec/nut          ire-gcrf
%
%  references    :
%    vallado       2004, 205-219
%
% [recef, vecef, aecef] = eci2ecef06(reci, veci, aeci, iau06arr, xysarr, ttt, jdut1, lod, xp, yp, ddx, ddy, opt1 )
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


