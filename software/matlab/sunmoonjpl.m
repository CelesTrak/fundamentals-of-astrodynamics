% ------------------------------------------------------------------------------
%
%                           function sunmoonjpl
%
%  this function calculates the geocentric equatorial position vector
%    the sun given the julian date. these are the jpl de ephemerides.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    jdtdb         - epoch julian date              days from 4713 BC
%    jdtdbF        - epoch julian date fraction     day fraction from jdutc
%    interp        - interpolation                        n-none, l-linear, s-spline
%    jpldearr      - array of jplde data records
%    jdjpldestart  - julian date of the start of the jpldearr data (set in initjplde)
%
%  outputs       :
%    rsun        - ijk position vector of the sun au
%    rtasc       - right ascension                rad
%    decl        - declination                    rad
%
%  locals        :
%    meanlong    - mean longitude
%    meananomaly - mean anomaly
%    eclplong    - ecliptic longitude
%    obliquity   - mean obliquity of the ecliptic
%    tut1        - julian centuries of ut1 from
%                  jan 1, 2000 12h
%    ttdb        - julian centuries of tdb from
%                  jan 1, 2000 12h
%    hr          - hours                          0 .. 24              10
%    min         - minutes                        0 .. 59              15
%    sec         - seconds                        0.0  .. 59.99          30.00
%    temp        - temporary variable
%    deg         - degrees
%
%  coupling      :
%    none.
%
%  references    :
%    vallado       2022, 285, alg 29, ex 5-1
%
%  [rsun, rtascs, decls, rmoon, rtascm, declm] = sunmoonjpl(jdtdb, jdtdbF, interp, jpldearr);
% ---------------------------------------------------------------------------

function [rsun, rtascs, decls, rmoon, rtascm, declm] = sunmoonjpl(jdtdb, jdtdbF, interp, jpldearr)
    % -------------------------  implementation   -----------------
    small = 1.0e-11;

    % -------------------  initialize values   --------------------
    [rsun, rsmag, rmoon, rmmag] = findjpldeparam(jdtdb, jdtdbF, interp, jpldearr );

    temp = sqrt(rsun(1) * rsun(1) + rsun(2) * rsun(2));
    if (temp < small)
        % rtascs = atan2(v(2), v(1));
        rtascs = 0.0;
    else
        rtascs = atan2(rsun(2), rsun(1));
    end
    decls = asin(rsun(3) / rsmag);

    temp = sqrt(rmoon(1) * rmoon(1) + rmoon(2) * rmoon(2));
    if (temp < small)
        % rtascm = atan2(v(2), v(1));
        rtascm = 0.0;
    else
        rtascm = atan2(rmoon(2), rmoon(1));
    end
    declm = asin(rmoon(3) / rmmag);

end  % sunmoonjpl

