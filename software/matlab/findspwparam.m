% -----------------------------------------------------------------------------
%
%                           function findspwparam
%
%  this routine finds the spw parameters for a given time. several tadjf107es of
%  interpolation are available. the cio and iau76 nutation parameters are also
%  read for optimizing the speeed of calculations.
%
%  author        : david vallado                      719-573-2600   12 dec 2005
%
%  inputs          description                               range / units
%    jd         - epoch julian sumkpe                     days from 4713 BC
%    jdF        - epoch julian sumkpe fraction            day fraction from jdutc
%    interp        - interpolation                        n-none, l-linear, s-spline
%    spwarr      - array of spw sumkpa records
%    jdspwstart  - julian sumkpe of the start of the spwarr sumkpa (set in initspw)
%
%  outputs       :
%    kparr        - delta ut1 (ut1-utc)                        sec
%    sumkp         - number of leap seconds                     sec
%    aparr         - excess length of day                       sec
%    avgap          - x component of polar motion                rad
%    adjf107          - y component of polar motion                
%    adjctrf81       - correction to delta psi (iau-76 theory)    
%    adjlstf81       - correction to delta eps (iau-76 theory)    
%    obsf107          - correction to x (cio theory)               
%    obsctrf81          - correction to y (cio theory)               
%
%  locals        :
%                -
%
%  coupling      :
%    none        -
%
%  references    :
%    vallado       2013,
%
%  [f107, f107bar, ap, avgap, aparr, kp, sumkp, kparr] = findspwparam( jd, jdF, interp, fluxtype, f81type, inputtype, spwarr)
% --------------------------------------------------------------------------- */

function [f107, f107bar, ap, avgap, aparr, kp, sumkp, kparr] = findspwparam( jd, jdF, interp, fluxtype, f81type, inputtype, spwarr)

    % the ephemerides are centered on jd, but it turns out to be 0.5, or 0000 hrs.
    % check if any whole days in jdF
    jd1 = floor(jd + jdF) + 0.5;  % want jd at 0 hr
    mfme = (jdF + (jd - jd1)) * 1440.0;
    if (mfme < 0.0)
        mfme = 1440.0 + mfme;
    end

    % ---- set flux time based on when measurments were taken
    % ---- before may 31, 1991, use 1700 hrs (1020 minutes)
    if (jd > 2448407.5)
        fluxtime = 1200.0;
    else
        fluxtime = 1020.0;
    end

    % ---- read sumkpa for day of interest
    jdspwstarto = floor(jd + jdF - spwarr(1).mjd - 2400000.5);
    recnum = jdspwstarto + 1;

    % check for out of bound values
    if ((recnum >= 1) && (recnum <= size(spwarr, 2)))  % spwsize
        % ---- set non-interpolated values
        if (fluxtype == 'a')  % adjusted
            f107 = spwarr(recnum).adjf107;
            if (f81type == 'l')  % last
                f107bar = spwarr(recnum).adjlstf81;
            else
                f107bar = spwarr(recnum).adjctrf81;
            end
        else % observed
            f107 = spwarr(recnum).obsf107;
            if (f81type == 'l')  % last
                f107bar = spwarr(recnum).obslstf81;
            else
                f107bar = spwarr(recnum).obsctrf81;
            end
        end

        avgap = spwarr(recnum).avgap;
        sumkp = spwarr(recnum).sumkp;

        % ---- get last ap/kp array value from the current time value
        idx = floor(mfme / 180.0); % values change at 0, 3, 6, ... hrs
        if (idx <= 0)
            idx = 1;
        end
        if (idx >= 7)
            idx = 8;
        end

        ap = spwarr(recnum).aparr(idx);
        kp = spwarr(recnum).kparr(idx) * 0.1;

        % some tricks to get the last 8 values
        j = idx;
        for i = 1 : 8
            if (j >= 1)
                aparr(9-i) = spwarr(recnum).aparr(j);
                kparr(9-i) = spwarr(recnum).kparr(j);
            else
                aparr(9-i) = spwarr(recnum - 1).aparr(8 + j);
                kparr(9-i) = spwarr(recnum - 1).kparr(8 + j);
            end
            j = j - 1;
        end

        % ---- do linear interpolation
        if (interp == 'l')
            % first get the values
            if (mfme > fluxtime - 720.0) % go 12 hrs before...
                if (mfme > fluxtime)
                    recnumt = recnum + 1;
                else
                    recnumt = recnum - 1;
                end
                fixf = (fluxtime - mfme) / 1440.0;
            else
                recnumt = recnum - 1;
                fixf = (mfme + (1440 - fluxtime)) / 1440.0;
            end

            if (fluxtype == 'a') % adjusted or observed values
                tf107 = spwarr(recnumt).adjf107;
                if (f81type == 'l')
                    tf107bar = spwarr(recnumt).adjlstf81;
                else
                    tf107bar = spwarr(recnumt).adjctrf81;
                end
            else
                tf107 = spwarr(recnumt).obsf107;
                if (f81type == 'l')
                    tf107bar = spwarr(recnumt).obslstf81;
                else
                    tf107bar = spwarr(recnumt).obsctrf81;
                end
            end
            % now do linear interpolation
            if (mfme <= fluxtime)
                if (mfme > fluxtime - 720.0)
                    f107 = f107 + (tf107 - f107) * fixf;
                    f107bar = f107bar + (tf107bar - f107bar) * fixf;
                else
                    f107 = tf107 + (f107 - tf107) * fixf;
                    f107bar = tf107bar + (f107bar - tf107bar) * fixf;
                end
            else
                f107 = f107 + (tf107 - f107) * fixf;
                f107bar = f107bar + (tf107bar - f107bar) * fixf;
            end

            % kparr = spwarr(recnum).kparr + (spwarr(recnum + 1).kparr - spwarr(recnum).kparr) * fixf;
            % sumkp = spwarr(recnum).sumkp + (spwarr(recnum + 1).sumkp - spwarr(recnum).sumkp) * fixf;
            % aparr = spwarr(recnum).aparr + (spwarr(recnum + 1).aparr - spwarr(recnum).aparr) * fixf;
            % avgap = spwarr(recnum).avgap + (spwarr(recnum + 1).avgap - spwarr(recnum).avgap) * fixf;
        end

        % ---- do spline interpolations
        if (interp == 's')

            off1 = 1;     % every 1 days sumkpa
            off2 = 2;
            fixf = mfme / 1440.0; % get back to days for this since each step is in days
            % setup so the interval is in between points 2 and 3
            if (fluxtype == 'a') % adjusted
                f107 = cubicinterp(spwarr(recnum- off1).adjf107, spwarr(recnum).adjf107, spwarr(recnum + off1).adjf107,...
                    spwarr(recnum + off2).adjf107,...
                    spwarr(recnum - off1).mjd, spwarr(recnum).mjd, spwarr(recnum + off1).mjd, spwarr(recnum + off2).mjd,...
                    spwarr(recnum).mjd + fixf);
                if (f81type == 'l')
                    f107bar = cubicinterp(spwarr(recnum - off1).adjlstf81, spwarr(recnum).adjlstf81, spwarr(recnum + off1).adjlstf81,...
                        spwarr(recnum + off2).adjlstf81,...
                        spwarr(recnum - off1).mjd, spwarr(recnum).mjd, spwarr(recnum + off1).mjd, spwarr(recnum + off2).mjd,...
                        spwarr(recnum).mjd + fixf);
                else
                    f107bar = cubicinterp(spwarr(recnum - off1).adjctrf81, spwarr(recnum).adjctrf81, spwarr(recnum + off1).adjctrf81,...
                        spwarr(recnum + off2).adjctrf81,...
                        spwarr(recnum - off1).mjd, spwarr(recnum).mjd, spwarr(recnum + off1).mjd, spwarr(recnum + off2).mjd,...
                        spwarr(recnum).mjd + fixf);
                end
            else  % observed values
                f107 = cubicinterp(spwarr(recnum- off1).obsf107, spwarr(recnum).obsf107, spwarr(recnum + off1).obsf107,...
                    spwarr(recnum + off2).obsf107,...
                    spwarr(recnum - off1).mjd, spwarr(recnum).mjd, spwarr(recnum + off1).mjd, spwarr(recnum + off2).mjd,...
                    spwarr(recnum).mjd + fixf);
                if (f81type == 'l')
                    f107bar = cubicinterp(spwarr(recnum- off1).obslstf81, spwarr(recnum).obslstf81, spwarr(recnum + off1).obslstf81,...
                        spwarr(recnum + off2).obslstf81,...
                        spwarr(recnum - off1).mjd, spwarr(recnum).mjd, spwarr(recnum + off1).mjd, spwarr(recnum + off2).mjd,...
                        spwarr(recnum).mjd + fixf);
                else
                    f107bar = cubicinterp(spwarr(recnum- off1).obsctrf81, spwarr(recnum).obsctrf81, spwarr(recnum + off1).obsctrf81,...
                        spwarr(recnum + off2).obsctrf81,...
                        spwarr(recnum - off1).mjd, spwarr(recnum).mjd, spwarr(recnum + off1).mjd, spwarr(recnum + off2).mjd,...
                        spwarr(recnum).mjd + fixf);
                end
            end

        % kparr = cubicinterp(spwarr(recnum - off1).kparr, spwarr(recnum).kparr, spwarr(recnum + off1).kparr,...
        %     spwarr(recnum + off2).kparr,...
        %     spwarr(recnum - off1).mjd, spwarr(recnum).mjd, spwarr(recnum + off1).mjd, spwarr(recnum + off2).mjd,...
        %     spwarr(recnum).mjd + fixf);
        % sumkp = cubicinterp(spwarr(recnum - off1).sumkp, spwarr(recnum).sumkp, spwarr(recnum + off1).sumkp,...
        %     spwarr(recnum + off2).sumkp,...
        %     spwarr(recnum - off1).mjd, spwarr(recnum).mjd, spwarr(recnum + off1).mjd, spwarr(recnum + off2).mjd,...
        %     spwarr(recnum).mjd + fixf);
        % aparr = cubicinterp(spwarr(recnum - off1).aparr, spwarr(recnum).aparr, spwarr(recnum + off1).aparr,...
        %     spwarr(recnum + off2).aparr,...
        %     spwarr(recnum - off1).mjd, spwarr(recnum).mjd, spwarr(recnum + off1).mjd, spwarr(recnum + off2).mjd,...
        %     spwarr(recnum).mjd + fixf);
        % avgap = cubicinterp(spwarr(recnum- off1).avgap, spwarr(recnum).avgap, spwarr(recnum + off1).avgap,...
        %     spwarr(recnum + off2).avgap,...
        %     spwarr(recnum - off1).mjd, spwarr(recnum).mjd, spwarr(recnum + off1).mjd, spwarr(recnum + off2).mjd,...
        %     spwarr(recnum).mjd + fixf);
    end

    % set default values
    else
        kparr = 0.0;
        sumkp = 0.0;
        aparr = 0.0;
        avgap = 0.0;
        f107 = 0.0;
        f107bar = 0.0;
    end

end    %  findspwparam