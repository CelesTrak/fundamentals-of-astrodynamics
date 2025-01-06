% -----------------------------------------------------------------------------
%
%                           function findeopparam
%
%  this routine finds the eop parameters for a given time. several types of
%  interpolation are available.
%
%  author        : david vallado                      719-573-2600   12 dec 2005
%
%  inputs          description                               range / units
%    jdtdb         - epoch julian date                     days from 4713 BC
%    jdtdbF        - epoch julian date fraction            day fraction from jdutc
%    interp        - interpolation                        n-none, l-linear, s-spline
%    eopearr      - array of eop data records
%    jdeopstart  - julian date of the start of the eopearr data (set in initeop)
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
%    vallado       2013,
%  [dut1, dat, lod, xp, yp, ddpsi, ddeps, ddx, ddy] = findeopparam( jdtdb, jdtdbF, interp, eoparr)
% --------------------------------------------------------------------------- 

function [dut1, dat, lod, xp, yp, ddpsi, ddeps, ddx, ddy] = findeopparam( jdtdb, jdtdbF, interp, eoparr)

    convrt = pi / (3600.0 * 180.0);  % " to rad
   
    % the ephemerides are centered on jdtdb, but it turns out to be 0.5, or 0000 hrs.
    % check if any whole days in jdF
    jdb = floor(jdtdb + jdtdbF) + 0.5;  % want jd at 0 hr
    mfme = (jdtdb + jdtdbF - jdb) * 1440.0;
    if (mfme < 0.0)
        mfme = 1440.0 + mfme;
    end

    % ---- read data for day of interest
    jdeopstarto = floor(jdtdb + jdtdbF - eoparr(1).mjd - 2400000.5);
    recnum = floor(jdeopstarto)+1;

    % check for out of bound values
    if ((recnum >= 1) && (recnum <= 51830))  % eopsize
        % ---- set non-interpolated values
        dut1 = eoparr(recnum).dut1;
        dat = eoparr(recnum).dat;
        lod = eoparr(recnum).lod;
        xp = eoparr(recnum).xp;
        yp = eoparr(recnum).yp;
        ddpsi = eoparr(recnum).ddpsi;
        ddeps = eoparr(recnum).ddeps;
        ddx = eoparr(recnum).dx;
        ddy = eoparr(recnum).dy;

        % ---- find nutation parameters for use in optimizing speed

        % ---- do linear interpolation
        if (interp == 'l')
            fixf = mfme / 1440.0;

            dut1 = eoparr(recnum).dut1 + (eoparr(recnum + 1).dut1 - eoparr(recnum).dut1) * fixf;
            dat = eoparr(recnum).dat + (eoparr(recnum + 1).dat - eoparr(recnum).dat) * fixf;
            lod = eoparr(recnum).lod + (eoparr(recnum + 1).lod - eoparr(recnum).lod) * fixf;
            xp = eoparr(recnum).xp + (eoparr(recnum + 1).xp - eoparr(recnum).xp) * fixf;
            yp = eoparr(recnum).yp + (eoparr(recnum + 1).yp - eoparr(recnum).yp) * fixf;
            ddpsi = eoparr(recnum).ddpsi + (eoparr(recnum + 1).ddpsi - eoparr(recnum).ddpsi) * fixf;
            ddeps = eoparr(recnum).ddeps + (eoparr(recnum + 1).ddeps - eoparr(recnum).ddeps) * fixf;
            ddx = eoparr(recnum).dx + (eoparr(recnum + 1).dx - eoparr(recnum).dx) * fixf;
            ddy = eoparr(recnum).dy + (eoparr(recnum + 1).dy - eoparr(recnum).dy) * fixf;
            %printf("sunm %i xp %lf fixf %lf n %lf nxt %lf \n", recnum, xp, fixf, eoparr(recnum).dut1, eoparr(recnum).dut1);
            %printf("recnum l %i fixf %lf %lf rsun %lf %lf %lf \n", recnum, fixf, eoparr(recnum).dut1, dut1, rsuny, rsunz);
        end

        % ---- do spline interpolations
        if (interp == 's')

            off1 = 1;     % every 1 days data
            off2 = 2;
            fixf = mfme / 1440.0; % get back to days for this since each step is in days
            % setup so the interval is in between points 2 and 3
            dut1 = cubicinterp(eoparr(recnum - off1).dut1, eoparr(recnum).dut1, eoparr(recnum + off1).dut1,...
                eoparr(recnum + off2).dut1,...
                eoparr(recnum - off1).mjd, eoparr(recnum).mjd, eoparr(recnum + off1).mjd, eoparr(recnum + off2).mjd,...
                eoparr(recnum).mjd + fixf);
            dat = cubicinterp(eoparr(recnum - off1).dat, eoparr(recnum).dat, eoparr(recnum + off1).dat,...
                eoparr(recnum + off2).dat,...
                eoparr(recnum - off1).mjd, eoparr(recnum).mjd, eoparr(recnum + off1).mjd, eoparr(recnum + off2).mjd,...
                eoparr(recnum).mjd + fixf);
            lod = cubicinterp(eoparr(recnum - off1).lod, eoparr(recnum).lod, eoparr(recnum + off1).lod,...
                eoparr(recnum + off2).lod,...
                eoparr(recnum - off1).mjd, eoparr(recnum).mjd, eoparr(recnum + off1).mjd, eoparr(recnum + off2).mjd,...
                eoparr(recnum).mjd + fixf);
            xp = cubicinterp(eoparr(recnum- off1).xp, eoparr(recnum).xp, eoparr(recnum + off1).xp,...
                eoparr(recnum + off2).xp,...
                eoparr(recnum - off1).mjd, eoparr(recnum).mjd, eoparr(recnum + off1).mjd, eoparr(recnum + off2).mjd,...
                eoparr(recnum).mjd + fixf);
            yp = cubicinterp(eoparr(recnum- off1).yp, eoparr(recnum).yp, eoparr(recnum + off1).yp,...
                eoparr(recnum + off2).yp,...
                eoparr(recnum - off1).mjd, eoparr(recnum).mjd, eoparr(recnum + off1).mjd, eoparr(recnum + off2).mjd,...
                eoparr(recnum).mjd + fixf);
            ddpsi = cubicinterp(eoparr(recnum - off1).ddpsi, eoparr(recnum).ddpsi, eoparr(recnum + off1).ddpsi,...
                eoparr(recnum + off2).ddpsi,...
                eoparr(recnum - off1).mjd, eoparr(recnum).mjd, eoparr(recnum + off1).mjd, eoparr(recnum + off2).mjd,...
                eoparr(recnum).mjd + fixf);
            ddeps = cubicinterp(eoparr(recnum - off1).ddeps, eoparr(recnum).ddeps, eoparr(recnum + off1).ddeps,...
                eoparr(recnum + off2).ddeps,...
                eoparr(recnum - off1).mjd, eoparr(recnum).mjd, eoparr(recnum + off1).mjd, eoparr(recnum + off2).mjd,...
                eoparr(recnum).mjd + fixf);
            ddx = cubicinterp(eoparr(recnum- off1).dx, eoparr(recnum).dx, eoparr(recnum + off1).dx,...
                eoparr(recnum + off2).dx,...
                eoparr(recnum - off1).mjd, eoparr(recnum).mjd, eoparr(recnum + off1).mjd, eoparr(recnum + off2).mjd,...
                eoparr(recnum).mjd + fixf);
            ddy = cubicinterp(eoparr(recnum- off1).dy, eoparr(recnum).dy, eoparr(recnum + off1).dy,...
                eoparr(recnum + off2).dy,...
                eoparr(recnum - off1).mjd, eoparr(recnum).mjd, eoparr(recnum + off1).mjd, eoparr(recnum + off2).mjd,...
                eoparr(recnum).mjd + fixf);
            %printf("recnum s %i mfme %lf days rsun %lf %lf %lf \n", recnum, mfme, rsunx, rsuny, rsunz);
            %printf(" %lf %lf %lf %lf \n", eoparr(recnum - off2).yp, eoparr(recnum - off1.yp, eoparr(recnum).yp, eoparr(recnum + off1).yp);
        end

        % set default values
    else
        dut1 = 0.0;
        dat = 0.0;
        lod = 0.0;
        xp = 0.0;
        yp = 0.0;
        ddpsi = 0.0;
        ddeps = 0.0;
        ddx = 0.0;
        ddy = 0.0;
    end

    % now convert units for use in operations
    xp = xp * convrt;  % " to rad
    yp = yp * convrt;
    ddpsi = ddpsi * convrt;  % " to rad
    ddeps = ddeps * convrt;
    ddx = ddx * convrt;  % " to rad
    ddy = ddy * convrt;

end    %  findeopparam