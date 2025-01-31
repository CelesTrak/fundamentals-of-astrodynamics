% -----------------------------------------------------------------------------
%
%                           function findxysparam
%
%  this routine finds the xys parameters for the iau 2006/2000 transformation.
%  several types of interpolation are available. this allows you to use the full cio
%  series, but maintain very fast performance.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    jdttt         - epoch julian date                      days from 4713 BC
%    jdftt         - epoch julian date fraction             day fraction from jdtt
%    interp        - interpolation                          n-none, l-linear, s-spline
%    xysarr        - array of xys data records              rad
%
%  outputs       :
%    x           - x component of cio                       rad
%    y           - y component of cio                       rad
%    s           -                                          rad
%
%  coupling      :
%    none        -
%
%  references    :
%
%  [x, y, s] = findxysparam( jdtt, jdttF, interp, xys06table);
% -----------------------------------------------------------------------------

function [x, y, s] = findxysparam( jdtt, jdttF, interp, xys06table)

    % check if any whole days in jdF
    jdb = floor(jdtt + jdttF) + 0.5;  % want jd at 0000 hr
    mfme = (jdtt + jdttF - jdb) * 1440.0;
    if (mfme < 0.0)
        mfme = 1440.0 + mfme;
    end

    % ---- read recnum for day of interest
    jdxysstarto = floor(jdtt + jdttF - xys06table.mjd_tt(1) - 2400000.5);
    recnum = floor(jdxysstarto) + 1;

    % check for out of bound values
    if ((recnum >= 1) && (recnum <= 51830))  % xyssize
        % ---- set non-interpolated values
        x = xys06table.x(recnum);
        y = xys06table.y(recnum);
        s = xys06table.s(recnum);

        % ---- do linear interpolation
        if (interp == 'l')
            fixf = mfme / 1440.0;

            x = xys06table.x(recnum) + (xys06table.x(recnum + 1) - xys06table.x(recnum)) * fixf;
            y = xys06table.y(recnum) + (xys06table.y(recnum + 1) - xys06table.y(recnum)) * fixf;
            s = xys06table.s(recnum) + (xys06table.s(recnum + 1) - xys06table.s(recnum)) * fixf;
        end

        % ---- do spline interpolations
        if (interp == 's')

            off1 = 1;     % every 1 days ya
            off2 = 2;
            fixf = mfme / 1440.0; % get back to days for this since each step is in days
            % setup so the interval is in between points 2 and 3
            x = cubicinterp(xys06table.x(recnum - off1), xys06table.x(recnum), xys06table.x(recnum + off1),...
                xys06table.x(recnum + off2),...
                xys06table.mjd_tt(recnum - off1), xys06table.mjd_tt(recnum), xys06table.mjd_tt(recnum + off1), xys06table.mjd_tt(recnum + off2),...
                xys06table.mjd_tt(recnum) + fixf);
            y = cubicinterp(xys06table.y(recnum - off1), xys06table.y(recnum), xys06table.y(recnum + off1),...
                xys06table.y(recnum + off2),...
                xys06table.mjd_tt(recnum - off1), xys06table.mjd_tt(recnum), xys06table.mjd_tt(recnum + off1), xys06table.mjd_tt(recnum + off2),...
                xys06table.mjd_tt(recnum) + fixf);
            s = cubicinterp(xys06table.s(recnum - off1), xys06table.s(recnum), xys06table.s(recnum + off1),...
                xys06table.s(recnum + off2),...
                xys06table.mjd_tt(recnum - off1), xys06table.mjd_tt(recnum), xys06table.mjd_tt(recnum + off1), xys06table.mjd_tt(recnum + off2),...
                xys06table.mjd_tt(recnum) + fixf);
        end

        % set default values
    else
        x = 0.0;
        y = 0.0;
        s = 0.0;
    end

end   %  findxysparam