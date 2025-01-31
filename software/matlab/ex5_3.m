    %     -----------------------------------------------------------------
    %
    %                              Ex5_3.m
    %
    %  this file demonstrates example 5-3.
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

    constmath;
    constastro;

    % --------  moon -------------- analytical moon ephemeris
    year = 1994;  % use UTC that will give TDT on the 28 Apr 0 hr tdt = utc + dat + 32.184
    mon = 4;
    day = 27;
    hr = 23;
    minute = 58;
    second = 59.816;
    dat = 28;
    dut1 = -0.0889721;

    [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac ] ...
        = convtime ( year, mon, day, hr, minute, second, timezone, dut1, dat );
    fprintf(1,'input data \n\n');
    fprintf(1,' year %5i ',year);
    fprintf(1,' mon %4i ',mon);
    fprintf(1,' day %3i ',day);
    fprintf(1,' %3i:%2i:%8.6f\n ',hr,minute,second );
    fprintf(1,' dut1 %8.6f s',dut1);
    fprintf(1,' dat %3i s',dat);
    fprintf(1,' ddpsi %8.6f " ddeps  %8.6f\n',ddpsi/conv, ddeps/conv);
    fprintf(1,'tt  %8.6f ttt  %16.12f jdtt  %18.11f \n',tt,ttt,jdtt + jdttfrac );
    fprintf(1,'tut1  %8.6f tut1  %16.12f jdut1  %18.11f \n',tut1,tut1,jdut1 + jdut1frac );

    [rmoon, rtasc, decl] = moon ( jd );
    
    fprintf(1,'rmoon  rtasc %14.6f deg decl %14.6f deg\n',rtasc*rad,decl*rad );
    fprintf(1,'rmoon new      %11.9f%11.9f%11.9f er\n',rmoon );
    fprintf(1,'rmoon new   %14.4f%14.4f%14.4f km\n',rmoon*re );

    %  1994  4 28  119377864.5535      84237799.9341      36522905.6699     150602138.5030  0.9867078     
    %  1994  4 28    -134038.3192       -311589.2121       -126061.1912
    rmoonaa = [-134038.3192  -311589.2121  -126061.1912 ];
    fprintf(1,'moon jplde430 -134038.3192  -311589.2121  -126061.1912 km\n' );

    
    
