% ----------------------------------------------------------------------------
%
%                           function tirs_eci
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
%    iau06arr    - iau2006 eop constants
%    ttt         - julian centuries of tt                        centuries
%    jdut1       - julian date of ut1                            days from 4713 bc
%    lod         - excess length of day                          sec
%    ddx         - delta x correction to gcrf                    rad
%    ddy         - delta y correction to gcrf                    rad
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
%    vallado       2022, 213
%
% [reci, veci, aeci] = tirs2eci(rtirs, vtirs, atirs, iau06arr, xysarr, ttt,
%            jdut1, lod, ddx, ddy, opt1 );
% ----------------------------------------------------------------------------

function [reci, veci, aeci] = tirs2eci(rtirs, vtirs, atirs, iau06arr, xysarr, ...
        ttt, jdut1, lod, ddx, ddy, opt1 )
    constastro;
    pnb = zeros(3,3);
    st = zeros(3,3);

    [fArgs06] = fundarg(ttt, '06');

    [prec, psia, wa, ea, xa] = precess ( ttt, '06' );

    % ---- ceo based, iau2006
    if not(contains(opt1, 'a')) || not(contains(opt1, 'b'))
        [x, y, s, pnb] = iau06xys (iau06arr, fArgs06, xysarr, ttt, ddx, ddy, opt1);
        [st, stdot] = sidereal(jdut1, 0.0, 0.0, 0.0, 0.0, 0, '06' );
    end

    % ---- class equinox based, 2000a
    if opt1 == '06a'
        [ deltapsi, pnb, prec, nut, l, l1, f, d, omega, ...
            lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate ...
            ] = iau06pna (ttt);
        [gst,st] = iau06gst(jdut1, ttt, deltapsi, l, l1, f, d, omega, ...
            lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate);
    end

    % ---- class equinox based, 2000b
    if opt1 == '06b'
        [ deltapsi, pnb, prec, nut, l, l1, f, d, omega, ...
            lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate ...
            ] = iau06pnb (ttt);
        [gst,st] = iau06gst(jdut1, ttt, deltapsi, l, l1, f, d, omega, ...
            lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate);
    end
    prec = eye(3);
    nut = pnb;


    thetasa= earthrot * (1.0  - lod/86400.0 );
    omegaearth = [0; 0; thetasa;];

    reci = prec*nut*st*rtirs;

    veci = prec*nut*st*(vtirs + cross(omegaearth,rtirs));

    temp = cross(omegaearth,rtirs);
    aeci = prec*nut*st*(atirs + cross(omegaearth,temp) ...
        + 2.0*cross(omegaearth,vtirs));

end