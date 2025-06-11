% ------------------------------------------------------------------------------
%
%                              Ex8_5.m
%
%  this file demonstrates example 8-5.
%
%                          companion code for
%             fundamentals of astrodyanmics and applications
%                                 2022
%                            by david vallado
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
% ------------------------------------------------------------------------------
    
    constastro;
    
        [normArr] = gravnorm(5+2);
    
    
    % LEO test
    recef = [-1033.4793830;  7901.2952754;  6380.3565958];
    vecef = [-3.225636520;  -2.872451450;   5.531924446];
    
    % GCRF newer version
    reci = [-605.79079600;	-5870.23042200;	3493.05191600];
    veci = [-1.568251000;	-3.702348000;	-6.479485000];
    aeci = [0.0; 0.0; 0.0];
    
    % j2000 values, original
    % reci = [-605.790430800; -5870.230407000; 3493.052004000];
    % veci = [-1.568251615; -3.702348353; -6.479484915];
    
    %reci = [487.0696937; -5330.5022406; 4505.7372146];
    %veci = [-2.101083975; 4.624581986; 5.688300377];
    
    % satellite charactereistics
    cd = 2.2;
    cr = 1.2;
    area = 40.0;  % m^2
    mass = 1000.0;  % kg
    
    
    conv = pi / (180.0*3600.0);
    
    % get time values
%    2020 02 18 58897  0.030745  0.336090 -0.1990821  0.0000677 -0.108022 -0.007466  0.000323 -0.000062  37

    year = 2020;
    mon = 2;
    day = 18;
    hr =  15;
    min=   8;
    sec=  47.23847;
    dut1 = -0.1991042;  % sec
    dat  = 37; % sec
    xp   =  0.030640 * conv;  % " to rad
    yp   =  0.336064 * conv;
    lod  =  0.0000643;
    ddpsi = -0.108193 * conv;  % " to rad
    ddeps = -0.007365 * conv;
    ddx = 0.000255 * conv;  % " to rad
    ddy = 0.000039 * conv;
    nutorder = 106;
    eqeterms = 2;
    timezone=0;
    
    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);
    
    normalized = 'y';
    fname = 'D:/Dataorig/Gravity/EGM2008 norm 486.txt';
    %fname = 'D:/Dataorig/Gravity/EGM2008_to2190_TideFree.txt';
    fprintf(1,'Reading %s \n', fname);
    [gravarr] = readgravityfield(fname, normalized);

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
    fprintf(1, ' order %3i  eqeterms %3i  \n', nutorder, eqeterms );
    fprintf(1, 'units are km and km/s and km/s2\n' );
    
    % -------- convtime    - convert time from utc to all the others
    fprintf(1, 'convtime results\n');
    [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] ...
        = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
    fprintf(1, 'ut1 %8.6f tut1 %16.12f jdut1 %18.11f ', ut1, tut1, jdut1+jdut1frac );
    [h, m, s] = sec2hms( ut1 );
    fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
    fprintf(1, 'utc %8.6f ', utc );
    [h, m, s] = sec2hms( utc );    fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
    fprintf(1, 'tai %8.6f', tai );
    [h, m, s] = sec2hms( tai );
    fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
    fprintf(1, 'tt  %8.6f ttt  %16.12f jdtt  %18.11f ', tt, ttt, jdtt+jdttfrac );
    [h, m, s] = sec2hms( tt );
    fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s);
    fprintf(1, 'tdb %8.6f ttdb %16.12f jdtdb %18.11f\n', tdb, ttdb, jdtdb+jdtdbfrac );
    
    
    % find ecef frame values for position and velocites
    
    [fArgs] = fundarg(ttt, '80');

    [prec,psia,wa,ea,xa] = precess ( ttt, '80' );

    [deltapsi, trueeps, meaneps, nut] = nutation  (ttt, ddpsi, ddeps, iau80arr, fArgs);
    
    [st, stdot] = sidereal(jdut1, deltapsi, meaneps, fArgs(5), lod, eqeterms, '80' );
    
    [pm] = polarm(xp,yp,ttt,'80');
    
    % gravity test
    recef = [ -2110.28953926559006504;  -5511.91602450908976607;  3491.91339855242995327];
    vecef = [ -2.8625582997;    -3.0199822163;    -6.4824698032]; 
    % two-body
    aecef = -mu*recef/(mag(recef)^3);  
    [reci, veci, aeci] = ecef2eci(recef, vecef, aecef, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );

    reci = [ -605.79120446;    -5870.23016547;    3493.05227627];
    veci = [-1.568252366;    -3.702348317;   -6.479484538];
    aeci = [7.483962e-04;     7.252325e-03;    -4.327587e-03];
 
    
    % -------- polarm      - transformation matrix for polar motion
    fprintf(1, 'polar motion matrix versions -------------------\n');
    [pm] = polarm(xp, yp, ttt, '80');
    fprintf(1, '%20.15f %20.15f %20.15f\n', pm' );
    
    
    thetasa= earthrot * (1.0  - lod/86400.0 );
    omegaearth = [0; 0; thetasa;];
    
    rpef  = st'*nut'*prec'*reci;
    recef = pm'*rpef;
    
    vpef  = st'*nut'*prec'*veci - cross( omegaearth,rpef );
    vecef = pm'*vpef;
    
    fprintf(1, '\nrecef  %14.7f %14.7f %14.7f\n', recef );
    fprintf(1, 'vecef  %14.9f %14.9f %14.9f\n', vecef );
    
    [recef2w, vecef2w, aecef2w] = eci2ecef(reci, veci, aeci, iau80arr, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps );
    fprintf(1, 'eci-ecef 2, w corr         %14.7f %14.7f %14.7f\n', recef2w );
    fprintf(1, ' v %14.9f %14.9f %14.9f\n', vecef2w );
    
    % find gravity potential acceleration
    % first find legendre polynomials
    
    % --------------------find latgc and lon---------------------- }
    [latgc, latgd, lon, hellp] = ecef2ll(recef);
    fprintf(1, 'lat lon  %14.9f %14.9f %14.9f\n', latgc*rad, lon*rad, hellp );
    
    magr = mag(recef);
    accel2bdf = -mu*recef/magr^3;
    
    fprintf(1,'--------------------Find Legendre polynomials -------------\n');
    
    % fprintf(1,'\nunnormalized coefficents --------------- \n');
    % fprintf(1,'  4  0 %16.11e  %16.11e \n', gravarr.c(5, 1), gravarr.s(5, 1) );
    % fprintf(1,'  4  4 %16.11e  %16.11e \n', gravarr.c(5, 5), gravarr.s(5, 5) );
    % fprintf(1,' 21  1 %16.11e  %16.11e \n', gravarr.c(22, 2), gravarr.s(22, 2)  );
    fprintf(1,'\nnormalized coefficents --------------- \n');
    fprintf(1,'  2  0 %16.11e  %16.11e \n', gravarr.cNor(2, 1), gravarr.sNor(2, 1) );
    fprintf(1,'  4  0 %16.11e  %16.11e \n', gravarr.cNor(5, 1), gravarr.sNor(5, 1) );
    fprintf(1,'  4  4 %16.11e  %16.11e \n', gravarr.cNor(5, 5), gravarr.sNor(5, 5) );
    fprintf(1,' 21  1 %16.11e  %16.11e \n', gravarr.cNor(22, 2), gravarr.sNor(22, 2)  );
    
    % [legarr, Legarr1] = legPoly(latgc, order);
    degree = 6;
    order = 6;
    [legarrGU] = legpolyGTDS(latgc, degree, order);
    [legarrMU, legarrMN] = legpolyMont (latgc, normArr, degree, order);

    % this seems to show that the Montenbruck approach has smaller
    % (1e-15) diffs than the gtds approach (1e-12)
    %legarrMU-legarrGU
    %legarrMN-legarrGN
    
    fprintf(1,'\nLegendre polynomials --------------- \n');
    
    str1 = '';
    str2 = '';
    for L = 1+1: 6
        str0 = num2str(L-1,'%3i ');
        str0 = [str0, ' 0 '];
        str1 = '';
        str2 = '';
        for m = 0+1: L
            %fprintf(1,' %3i  %3i  %16.11e ', L, m, legarr(L, m));
            %fprintf(1,' %3i  %3i  %16.11e ', L, m, legarr1(L, m));
            %fprintf(1,' %3i  %3i  %16.11e ', L, m, legarrEx(L, m));
            %fprintf(1,' %3i  %3i  %16.11e ', L, m, legarr(L, m));
            s = num2str(legarrMN(L, m),'%19.10e ');
            str1 = [str1, ' ', s];
        end
        fprintf(1,' %s %s \n', str0, str1);
        fprintf(1,' %s %s \n\n', str0, str2);
    end
    
    [trigarr] = trigpolyGTDS(latgc, lon, degree+3);
    % [VArr, WArr] = trigpolyMont(recef, latgc, degree+3);
   
    % [aPertGTDS] = GravAccelGTDS ( recef, gravarr, normArr, degree, order);

    [LegGottN, accel] = GravAccelGott(recef, gravarr, degree, order);
    fprintf(1,'gott R ecef wo2b     %2d %2d %21.17f  %21.17f  %21.17f  %21.17f \n', ...
        degree, order, accel' - accel2bdf, mag(accel' - accel2bdf));
    fprintf(1, ' accel w2b  %14.9f %14.9f %14.9f\n\n', accel );
    
    
    
    fprintf(1,'------------------ find drag acceleration\n');
    density = 1.5e-12 ;  % kg/m3
    magv = mag(vecef);
    vrel(1) = vecef(1);  % vecef normal is veci to tod, then - wxr
    vrel(2) = vecef(2);
    vrel(3) = vecef(3);
    fprintf(1, ' vrel %14.9f %14.9f %14.9f\n', vrel );
    adrag = -0.5 * density * cd * area/mass * magv * vrel * 1000.0 ;  % simplify vel, get units to km
    
    fprintf(1, ' adrag efc %16.11f %16.11f %16.11f\n\n', adrag );
    aeci = prec*nut*st*pm*adrag';
    fprintf(1, ' adrag eci %16.11f %16.11f %16.11f\n', aeci );
    
    
    fprintf(1,' ------------------ find third body acceleration\n');
    % sun
    infilename = append('D:\Codes\LIBRARY\DataLib\', 'sunmooneph_430t12.txt');
    [jpldearr] = readjplde(infilename);

    [rsun, rsmag, rmoon, rmmag] = findjpldeparam(jdtdb+1, jdtdbfrac, 's', jpldearr);
    % hard code to stk
    rsuns = [126916355.384390;    -69567131.339884;    -30163629.424510]';
    % JPL ans  2020  2 18  M          0.6306
    rmoonj = [14462.2967; -357096.9762; -151599.3021]';
    %JPL ans  2020  2 18 15:08:47.23847 S       0.6306
    rsunj = [126921698.4134; -69564121.8695; -30156263.9220]';
    
    fprintf(1, ' diff rsun stk-mine %14.5f %14.5f %14.5f %14.5f\n', rsuns-rsun, mag(rsuns-rsun) );
    fprintf(1, ' diff rsun jpl-mine %14.5f %14.5f %14.5f %14.5f\n', rsunj-rsun, mag(rsunj-rsun) );
    fprintf(1, ' diff rsun stk-jpl  %14.5f %14.5f %14.5f %14.5f\n', rsuns-rsunj,mag(rsuns-rsunj) );
    fprintf(1, ' diff rmoon jpl-mine %14.5f %14.5f %14.5f %14.5f\n', rmoonj-rmoon,mag(rmoonj-rmoon) );
    fprintf(1, ' rsun  %14.5f %14.5f %14.5f\n', rsun );
    fprintf(1, ' rmoon %14.9f %14.9f %14.9f\n', rmoon );
    
    mu3 = musun;
    rsat3 = rsun - reci';
    magrsat3 = mag(rsat3);
    rearth3 = rsun;
    magrearth3 = mag(rearth3);
    athirdbody = mu3*(rsat3/magrsat3^3 - rearth3/magrearth3^3);
    fprintf(1, ' a3bodyS  eci %16.11f %16.11f %16.11f\n', athirdbody );
    athirdbody2 = -mu3/magrearth3^3 * (rearth3 - 3.0*rearth3 * (dot(reci', rearth3)/magrearth3^2) ...
        - 7.5 * rearth3 * ((dot(reci', rearth3)/magrearth3^2))^2 );
    fprintf(1, ' a3bodyS2 eci %16.11f %16.11f %16.11f\n', athirdbody2 );
    q = (mag(reci)^2 + 2.0 * dot(reci', rsat3)) * (magrearth3^2 + magrearth3*magrsat3 + magrsat3^2) / (magrearth3^3*magrsat3^3*(magrearth3 + magrsat3));
    athirdbody1 =  mu3*(rsat3*q - reci'/magrearth3^3);
    fprintf(1, ' a3bodyS1 eci %16.11f %16.11f %16.11f\n', athirdbody1 );
    a3body = athirdbody1;
    
    % moon
    mu3 = mumoon;
    rsat3 = rmoon - reci';
    magrsat3 = mag(rsat3);
    rearth3 = rmoon;
    magrearth3 = mag(rearth3);
    athirdbody = mu3*(rsat3/magrsat3^3 - rearth3/magrearth3^3);
    fprintf(1, ' a3bodyM  eci %16.11f %16.11f %16.11f\n\n', athirdbody );
    athirdbody2 = -mu3/magrearth3^3 * (rearth3 - 3.0*rearth3 * (dot(reci', rearth3)/magrearth3^2) ...
        - 7.5 * rearth3 * (dot(reci', rearth3)/magrearth3^2) );
    fprintf(1, ' a3bodyM2 eci %16.11f %16.11f %16.11f\n', athirdbody2 );
    q = (mag(reci)^2 + 2.0 * dot(reci', rsat3)) * (magrearth3^2 + magrearth3*magrsat3 + magrsat3^2) / (magrearth3^3*magrsat3^3*(magrearth3 + magrsat3));
    athirdbody1 =  mu3*(rsat3*q - reci'/magrearth3^3);
    fprintf(1, ' a3bodyM1 eci %16.11f %16.11f %16.11f\n', athirdbody1 );
    a3body = a3body + athirdbody1;
    
    fprintf(1,' ------------------ find srp acceleration\n');
    psrp = 4.57e-6;  % N/m^2 = kgm/s2 / m2
    rsatsun = rsun - reci';
    magrsatsun = mag(rsatsun);
    asrp = -(psrp * cr * area/mass * rsatsun/magrsatsun) / 1000.0;
    fprintf(1, ' asrp eci %16.11f %16.11f %16.11f\n\n', asrp );
    
    fprintf(1,' add perturbing accelerations\n');
    % add 3body for both?
    aecef =  adrag + athirdbody + asrp;  % agravity +
    fprintf(1, ' aecef %16.11f %16.11f %16.11f\n', aecef );
    
    % ---- move acceleration from earth fixed coordinates to eci
    % there are no cross products here as normal
    aeci = prec*nut*st*pm*aecef';
    fprintf(1, ' aeci %16.11f %16.11f %16.11f\n', aeci );
    
    % find two body component of eci acceleration
    aeci2 = -mu*reci/(mag(reci)^3);
    fprintf(1, ' aeci2 %16.11f %16.11f %16.11f\n', aeci2 );
    
    % total acceleration
    aeci = aeci2 + aeci;
    fprintf(1, 'total aeci %16.11f %16.11f %16.11f\n', aeci );
    
    
    
    

