% -----------------------------------------------------------------------------
%
%                           function findjpldeparam
%
%  this routine finds the jplde parameters for a given time. several types of
%  interpolation are available. the cio and iau76 nutation parameters are also
%  read for optimizing the speeed of calculations.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    jdtdb         - epoch julian date                      days from 4713 BC
%    jdtdbF        - epoch julian date fraction             day fraction from jdutc
%    interp        - interpolation                          n-none, l-linear, s-spline
%    jpldearr      - array of jplde data records
%    jdjpldestart  - julian date of the start of the jpldearr data (set in initjplde)
%
%  outputs       :
%    dut1        - delta ut1 (ut1-utc)                        sec
%    dat         - number of leap seconds                     sec
%    lod         - excess length of day                       sec
%    xp          - x component of polar motion                rad
%    yp          - y component of polar motion                rad
%    ddpsi       - correction to delta psi (iau-76 theory)    rad
%    ddeps       - correction to delta eps (iau-76 theory)    rad
%    dx          - correction to x (cio theory)               rad
%    dy          - correction to y (cio theory)               rad
%    x           - x component of cio                         rad
%    y           - y component of cio                         rad
%    s           -                                            rad
%    deltapsi    - nutation longitude angle                   rad
%    deltaeps    - obliquity of the ecliptic correction       rad
%
%  locals        :
%                -
%
%  coupling      :
%    none        -
%
%  references    :
%    vallado       2022, 612
%
%  [rsun, rsmag, rmoon, rmmag] = findjpldeparam( jdtdb, jdtdbF, 'l', jpldearr);
% ---------------------------------------------------------------------------

function [rsun, rsmag, rmoon, rmmag] = findjpldeparam( jdtdb, jdtdbF, interp, jpldearr)
    % the ephemerides are centered on jdtdb, but it turns out to be 0.5, or 0000 hrs.
    % check if any whole days in jdF
    jdb = floor(jdtdb + jdtdbF) + 0.5;  % want jd at 0 hr
    mfme = (jdtdb + jdtdbF - jdb) * 1440.0;
    if (mfme < 0.0)
        mfme = 1440.0 + mfme;
    end

    %printf("jdtdb %lf  %lf  %lf  %lf \n ", jdtdb, jdtdbF, jdb, mfme);x[recnum]

    % ---- read data for day of interest
    jdjpldestarto = floor(jdtdb + jdtdbF - jpldearr.mjd(1) - 2400000.5);
    recnum = floor(jdjpldestarto) * 2 + 3;

    % check for out of bound values
    if ((recnum >= 1) && (recnum <= 51830))  % jpldesize
        % ---- set non-interpolated values
        rsun(1) = jpldearr.rsun1(recnum);
        rsun(2) = jpldearr.rsun2(recnum);
        rsun(3) = jpldearr.rsun3(recnum);
        rsmag = jpldearr.rsmag;
        rmoon(1) = jpldearr.rmoon1(recnum);
        rmoon(2) = jpldearr.rmoon2(recnum);
        rmoon(3) = jpldearr.rmoon3(recnum);
        rmmag = mag(rmoon);

        % ---- find nutation parameters for use in optimizing speed

        % ---- do linear interpolation
        if (interp == 'l')
            fixf = mfme / 1440.0;

            rsun(1) = jpldearr.rsun1(recnum) + (jpldearr.rsun1(recnum + 1) - jpldearr.rsun1(recnum)) * fixf;
            rsun(2) = jpldearr.rsun2(recnum) + (jpldearr.rsun2(recnum + 1) - jpldearr.rsun2(recnum)) * fixf;
            rsun(3) = jpldearr.rsun3(recnum) + (jpldearr.rsun3(recnum + 1) - jpldearr.rsun3(recnum)) * fixf;
            rsmag = mag(rsun);
            rmoon(1) = jpldearr.rmoon1(recnum) + (jpldearr.rmoon1(recnum + 1) - jpldearr.rmoon1(recnum)) * fixf;
            rmoon(2) = jpldearr.rmoon2(recnum) + (jpldearr.rmoon2(recnum + 1) - jpldearr.rmoon2(recnum)) * fixf;
            rmoon(3) = jpldearr.rmoon3(recnum) + (jpldearr.rmoon3(recnum + 1) - jpldearr.rmoon3(recnum)) * fixf;
            rmmag = mag(rmoon);
            %printf("sunm %i rsmag %lf fixf %lf n %lf nxt %lf \n", recnum, rsmag, fixf, jpldearr.rsun(1), jpldearr.rsun(1));
            %printf("recnum l %i fixf %lf %lf rsun %lf %lf %lf \n", recnum, fixf, jpldearr.rsun(1), rsun(1), rsuny, rsunz);
        end

        % ---- do spline interpolations
        if (interp == 's')
            off1 = 1;     % every 1 days data
            off2 = 2;
            fixf = mfme / 1440.0; % get back to days for this since each step is in days
            % setup so the interval is in between points 2 and 3
            rsun(1) = cubicinterp(jpldearr.rsun1(recnum - off1), jpldearr.rsun1(recnum), jpldearr.rsun1(recnum + off1),...
                jpldearr.rsun1(recnum + off2),...
                jpldearr.mjd(recnum - off1), jpldearr.mjd(recnum), jpldearr.mjd(recnum + off1), jpldearr.mjd(recnum + off2),...
                jpldearr.mjd(recnum) + fixf);
            rsun(2) = cubicinterp(jpldearr.rsun2(recnum - off1), jpldearr.rsun2(recnum), jpldearr.rsun2(recnum + off1),...
                jpldearr.rsun2(recnum + off2),...
                jpldearr.mjd(recnum - off1), jpldearr.mjd(recnum), jpldearr.mjd(recnum + off1), jpldearr.mjd(recnum + off2),...
                jpldearr.mjd(recnum) + fixf);
            rsun(3) = cubicinterp(jpldearr.rsun3(recnum - off1), jpldearr.rsun3(recnum), jpldearr.rsun3(recnum + off1),...
                jpldearr.rsun3(recnum + off2),...
                jpldearr.mjd(recnum - off1), jpldearr.mjd(recnum), jpldearr.mjd(recnum + off1), jpldearr.mjd(recnum + off2),...
                jpldearr.mjd(recnum) + fixf);
            rsmag = mag(rsun);
            rmoon(1) = cubicinterp(jpldearr.rmoon1(recnum - off1), jpldearr.rmoon1(recnum), jpldearr.rmoon1(recnum + off1),...
                jpldearr.rmoon1(recnum + off2),...
                jpldearr.mjd(recnum - off1), jpldearr.mjd(recnum), jpldearr.mjd(recnum + off1), jpldearr.mjd(recnum + off2),...
                jpldearr.mjd(recnum) + fixf);
            rmoon(2) = cubicinterp(jpldearr.rmoon2(recnum - off1), jpldearr.rmoon2(recnum), jpldearr.rmoon2(recnum + off1),...
                jpldearr.rmoon2(recnum + off2),...
                jpldearr.mjd(recnum - off1), jpldearr.mjd(recnum), jpldearr.mjd(recnum + off1), jpldearr.mjd(recnum + off2),...
                jpldearr.mjd(recnum) + fixf);
            rmoon(3) = cubicinterp(jpldearr.rmoon3(recnum- off1), jpldearr.rmoon3(recnum), jpldearr.rmoon3(recnum + off1),...
                jpldearr.rmoon3(recnum + off2),...
                jpldearr.mjd(recnum - off1), jpldearr.mjd(recnum), jpldearr.mjd(recnum + off1), jpldearr.mjd(recnum + off2),...
                jpldearr.mjd(recnum) + fixf);
            rmmag = mag(rmoon);
            %printf("recnum s %i mfme %lf days rsun %lf %lf %lf \n", recnum, mfme, rsunx, rsuny, rsunz);
            %printf(" %lf %lf %lf %lf \n", jpldearr(recnum - off2).mjd, jpldearr(recnum - off1.mjd, jpldearr(recnum).mjd, jpldearr(recnum + off1).mjd);
        end

        % set default values
    else
        rsun(1) = 0.0;
        rsun(2) = 0.0;
        rsun(3) = 0.0;
        rsmag = 0.0;
        rmoon(1) = 0.0;
        rmoon(2) = 0.0;
        rmoon(3) = 0.0;
        rmmag = 0.0;
    end

end    %  findjpldeparam