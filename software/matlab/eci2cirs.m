% ----------------------------------------------------------------------------
%
%                           function eci_cirs
%
%  this function transforms between the eci mean equator mean equinox (gcrf), and
%    the true of date frame (tod).
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    reci        - position vector eci                           km
%    veci        - velocity vector eci                           km/s
%    enum        - direction                                     eto, efrom
%    iau06arr    - iau2006 eop constants
%    jdtt        - julian date of tt                             days from 4713 bc
%    jdftt       - fractional julian date of tt                  days
%    ddx         - delta x correction to eci                     rad
%    ddy         - delta y correction to eci                     rad
%
%  outputs       :
%    rcirs       - position vector cirs                          km
%    vcirs       - velocity vector cirs                          km/s
%
%  locals        :
%    deltapsi    - nutation angle                                rad
%    trueeps     - true obliquity of the ecliptic                rad
%    meaneps     - mean obliquity of the ecliptic                rad
%    prec        - matrix for mod - eci
%    nut         - matrix for tod - mod
%
%  coupling      :
%   precess      - rotation for precession
%   nutation     - rotation for nutation
%
%  references    :
%    vallado       2022, 223-231
%
% [rcirs, vcirs, acirs, pnb] = eci2cirs( reci, veci, aeci, iau06arr, ttt, ddx, ddy, opt1 );
% ----------------------------------------------------------------------------

function [rcirs, vcirs, acirs] = eci2cirs( reci, veci, aeci, iau06arr, xysarr, ttt, ddx, ddy, opt1 )

    %      sethelp;

    [fArgs06] = fundarg(ttt, '06');

    % ---- ceo based, iau2000
    if not(contains(opt1, 'a')) || not(contains(opt1, 'b'))
        [x, y, s, pnb] = iau06xys (iau06arr, fArgs06, xysarr, ttt, ddx, ddy, opt1);
    end

    % ---- class equinox based, 2000a
    if opt1 == '06a'
        [ deltapsi, pnb, prec, nut, l, l1, f, d, omega, ...
            lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate ...
            ] = iau06pna (ttt);
    end

    % ---- class equinox based, 2000b
    if opt1 == '06b'
        [ deltapsi, pnb, prec, nut, l, l1, f, d, omega, ...
            lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate ...
            ] = iau06pnb (ttt);
    end

    rcirs = pnb'*reci;
    vcirs = pnb'*veci;
    acirs = pnb'*aeci;

end