% ----------------------------------------------------------------------------
%
%                           function cirs_ecef
%
%  this function transforms a vector from earth fixed (itrf) frame to
%  true of date (tod).
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    recef       - position vector earth fixed                   km
%    vecef       - velocity vector earth fixed                   km/s
%    direct      - direction of transfer                         eto, efrom
%    iau06arr    - iau2006 eop constants
%    ttt         - julian centuries of tt                        centuries
%    jdut1       - julian date of ut1                            days from 4713 bc
%    lod         - excess length of day                          sec
%    xp          - polar motion coefficient                      rad
%    yp          - polar motion coefficient                      rad
%    ddx         - delta x correction to gcrf                    rad
%    ddy         - delta y correction to gcrf                    rad
%
%  outputs       :
%    rcirs       - position vector eci                           km
%    vcirs       - velocity vector eci                           km/s
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
%    vallado       2022, 227
%
% [recef,vecef,aecef] = cirs2ecef(rcirs, vcirs, acirs, iau06arr, ttt, jdut1, ...
%       lod, xp, yp, ddx, ddy, opt1 );
% ----------------------------------------------------------------------------

function [recef,vecef,aecef] = cirs2ecef(rcirs, vcirs, acirs, iau06arr, ttt, ...
        jdut1, lod, xp, yp, ddx, ddy, opt1 )
    constastro;
    %      sethelp;

    [fArgs06] = fundarg(ttt, '06');

    % ---- ceo based, iau2000
    if contains(opt1, '06')
        [st, stdot] = sidereal(jdut1, 0.0, 0.0, 0.0, 0.0, 0, '06');
    end

    % ---- class equinox based, 2000a
    if opt1 == '06a'
        [gst,st] = iau00gst(jdut1, ttt, deltapsi, l, l1, f, d, omega, ...
            lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate);
    end

    % ---- class equinox based, 2000b
    if opt1 == '06b'
        [gst,st] = iau00gst(jdut1, ttt, deltapsi, l, l1, f, d, omega, ...
            lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate);
    end

    [pm] = polarm(xp,yp,ttt,'06');

    % ---- setup parameters for velocity transformations
    thetasa= earthrot * (1.0  - lod/86400.0 );
    omegaearth = [0; 0; thetasa;];

    rpef  = st'*rcirs;
    recef = pm'*rpef;

    vpef  = st'*vcirs - cross( omegaearth,rpef );
    vecef = pm'*vpef;

    temp  = cross(omegaearth,rpef);
    aecef = pm'*(st'*acirs - cross(omegaearth,temp) - 2.0*cross(omegaearth,vpef));

end