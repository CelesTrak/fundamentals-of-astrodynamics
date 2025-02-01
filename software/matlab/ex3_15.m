% ------------------------------------------------------------------------------
%
%                              Ex3_15.m
%
%  this file demonstrates example 3-15.
%
%                          companion code for
%             fundamentals of astrodyanmics and applications
%                                 2022
%                            by david vallado
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
% ------------------------------------------------------------------------------


    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

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
    fprintf(1, ' order %3i  eqeterms %3i  \n', order, eqeterms );
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

    % J2000
    [recii, vecii, aecii] = ecef2eci(recef, vecef, aecef, iau80arr, ttt, jdut1+jdut1frac, lod, xp, yp, eqeterms, 0.0, 0.0 );
    fprintf(1,'J2000 wo corr IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recii(1), recii(2), recii(3), vecii(1), vecii(2), vecii(3));

    % GCRF
    [reci, veci, aeci] = ecef2eci(recef, vecef, aecef, iau80arr, ttt, jdut1+jdut1frac, lod, xp, yp, eqeterms, ddpsi, ddeps );
    fprintf(1,'GCRF w corr   IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', reci(1), reci(2), reci(3), veci(1), veci(2), veci(3));
  
    [recefi, vecefi, aecefi] = eci2ecef(reci, veci, aeci, iau80arr, ttt, jdut1+jdut1frac, lod, xp, yp, eqeterms, ddpsi, ddeps );
    fprintf(1,'ITRF rev      IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', recefi(1), recefi(2), recefi(3), vecefi(1), vecefi(2), vecefi(3));

    