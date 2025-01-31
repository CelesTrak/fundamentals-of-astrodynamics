    %     -----------------------------------------------------------------
    %
    %                              ex3_14.m
    %
    %  this file demonstrates example 3-14.
    %
    %                          companion code for
    %             fundamentals of astrodynamics and applications
    %                                 2013
    %                            by david vallado
    %
    %     (h)               email davallado@gmail.com
    %     (w) 719-573-2600, email dvallado@agi.com
    %
    %     *****************************************************************
    %
    %  current :
    %            16 feb 19  david vallado
    %                         update for new constants
    %  changes :
    %            13 feb 07  david vallado
    %                         original baseline
    %
    %     *****************************************************************

    % LEO test
    recef = [-1033.4793830;  7901.2952754;  6380.3565958];
    vecef = [-3.225636520;  -2.872451450;   5.531924446];
    aecef = [0.001;0.002;0.003];

    conv = pi / (180.0*3600.0);

    year=2004;
    mon = 4;
    day = 6;
    hr =  7;
    min= 51;
    sec= 28.386009;

    dut1 = -0.4399619;  % sec
    dat  = 32;         % sec
    xp   = -0.140682 * conv;  % " to rad
    yp   =  0.333309 * conv;
    lod  =  0.0015563;
    ddpsi = -0.052195 * conv;  % " to rad
    ddeps = -0.003875 * conv;
    ddx = -0.000205 * conv;  % " to rad
    ddy = -0.000136 * conv;
    order = 106;
    eqeterms = 2;
    timezone=0;
    opt = 'c'; % specify the iau00 cio approach

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau06arr] = iau06in(fileLoc);

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [xys06table] = readxys(fileLoc);


    fprintf(1, 'input data \n\n');
    fprintf(1, ' year %5i ', year);
    fprintf(1, ' mon %4i ', mon);
    fprintf(1, ' day %3i ', day);
    fprintf(1, ' %3i:%2i:%8.6f\n', hr, min, sec );
    fprintf(1, ' dut1 %8.6f s', dut1);
    fprintf(1, ' dat %3i s', dat);
    fprintf(1, ' xp %8.6f "', xp / conv);
    fprintf(1, ' yp %8.6f "', yp / conv);
    fprintf(1, ' lod %8.6f s\n', lod);
    fprintf(1, ' ddpsi %8.6f " ddeps  %8.6f\n', ddpsi/conv, ddeps/conv);
    fprintf(1, ' ddx   %8.6f " ddy    %8.6f\n', ddx/conv, ddy/conv);
    fprintf(1, ' order %3i  eqeterms %3i  opt %3s \n', order, eqeterms, opt );
    fprintf(1, 'units are km and km/s and km/s2\n' );

    % -------- convtime    - convert time from utc to all the others
    fprintf(1, 'convtime results\n');
    [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] ...
        = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
    fprintf(1, 'ut1 %8.6f tut1 %16.12f jdut1 %18.11f ', ut1, tut1, jdut1+jdut1frac );
    [h, m, s] = sec2hms( ut1 );
    fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
    fprintf(1, 'utc %8.6f ', utc );
    [h, m, s] = sec2hms( utc );
    fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
    fprintf(1, 'tai %8.6f', tai );
    [h, m, s] = sec2hms( tai );
    fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
    fprintf(1, 'tt  %8.6f ttt  %16.12f jdtt  %18.11f ', tt, ttt, jdtt+jdttfrac );
    [h, m, s] = sec2hms( tt );
    fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
    fprintf(1, 'tdb %8.6f ttdb %16.12f jdtdb %18.11f\n', tdb, ttdb, jdtdb+jdtdbfrac );

    [recigg vecigg, aecig] = ecef2eci06(recef, vecef, aecef, iau06arr, xys06table, ...
        ttt, jdut1+jdut1frac, lod, xp, yp, ddx, ddy, '06x');
    fprintf(1, 'GCRF  IAU-2006 CIO %14.7f %14.7f %14.7f', recigg );
    fprintf(1, ' v %14.9f %14.9f %14.9f', vecigg );
    fprintf(1, ' a %14.9f %14.9f %14.9f\n', aecig );

