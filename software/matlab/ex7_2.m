    %     -----------------------------------------------------------------
    %
    %                              Ex7_2.m
    %
    %  this file demonstrates example 7-2.
    %
    %                          companion code for
    %             fundamentals of astrodynamics and applications
    %                                 2007
    %                            by david vallado
    %
    %     (w) 719-573-2600, email dvallado@agi.com
    %
    %     *****************************************************************
    %
    %  current :
    %             9 oct 07  david vallado
    %                         original
    %  changes :
    %             9 oct 07  david vallado
    %                         original baseline
    %
    %     *****************************************************************
    
    constmath;
    constastro;
    
    %     re = 6378.137;
    %     mu = 3.986004418e5;
    tu = 86400.0;
    
    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    %    re = 6378.145;
    %    mu = 3.986005e5;
    
    %    re = 149597870.0;  % km in 1 au
    %    mu = 1.32712428e11;
    %    tu = 86400.0;
    
    %    re = 1.0;  % 1 au
    %    mu = 1.0;
    %    tu = 1.0 / 58.132440906; % days in one solar tu
    
    convrt = pi / (180*3600.0);
    
    casenum = 0;  % 0 book
    %  casenum = 4;  % test
    
    diffsites = 'n';
    
    % typerun = 'l'; % laplace
    % typerun = 'g'; % gauss
    typerun = 'd'; % doubler
    % typerun = 'o'; % gooding
    %typerun = 'a';  % all
    
    %        filedat =load(append('Sat11Ex',int2str(casenum),'.dat'));
    
    
    % loop through all the test cases in the file
    for casenum = 0:30
        % ------ read the file line by line since a mix of text and data
        fid = fopen('d:\codes\library\datalib\anglestest.dat');
        tline = fgetl(fid);
        filedat = zeros(15, 11);
        ktr = 1;
        numcase = 0;
        numobs = 0;
        while ischar(tline)
            % check for a case number
            if (tline(1) == '#')
                %casenum = str2num(tline(2:4));
                answer = tline;
                ktr = 1;
                numcase = numcase + 1;
                oldnumobs = numobs;
                numobs = 0;
            else
                strn = regexprep(tline,',',' ');
                vars = sscanf(strn,'%d %d %d %d %d %f %f %f %f %f %f %s'); %# sscanf can read only numeric data :(
                filedat(ktr, 3) = vars(3);
                filedat(ktr, 2) = vars(2);
                filedat(ktr, 1) = vars(1);
                filedat(ktr, 4) = vars(4);
                filedat(ktr, 5) = vars(5);
                filedat(ktr, 6) = vars(6);
                filedat(ktr, 7) = vars(7);  % lat
                filedat(ktr, 8) = vars(8);  % lon
                filedat(ktr, 9) = vars(9);  % alt km
                filedat(ktr, 10) = vars(10); % rtasc
                filedat(ktr, 11) = vars(11); % decl
    
                numobs = numobs + 1;
                ktr = ktr + 1;
            end
            tline = fgetl(fid);
            if (~ischar(tline) || numcase == casenum+3)
                break;
            end   % end of file
        end
        fclose(fid);
        % numcase
        numobs = oldnumobs;
    
        obs1 = 1;
        obs2 = 2;
        obs3 = 3;
    
    
        %      mfme = currobsrec.hr(j) * 60 + currobsrec.min(j) + currobsrec.sec(j)/60.0;
        %      findeopparam ( currobsrec.jd, mfme, interp, eoparr, jdeopstart,
        %                     dut1, dat, lod, xp, yp, ddpsi, ddeps,
        %                     iaudx, dy, icrsx, y, s, deltapsi, deltaeps );
        %      [ut1, tut1, jdut1, utc, tai, tt, ttt, jdtt, tdb, ttdb, jdtdb ] ...
        %               = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
    
    
        if casenum == 0
            r2ans = [5897.954130507     5791.046114526     6682.733686585];
            v2ans = [  -4.393910234     4.576816355     1.482423676];
    
            dut1 = -0.609641;
            dat   = 35;
            xp    =  0.137495 * convrt;
            yp    =  0.342416 * convrt;
            lod   =  0.0;
            timezone= 0;
            terms = 0;
            ddpsi = 0.0;
            ddeps = 0.0;
            obs1 = 3;
            obs2 = 5;
            obs3 = 6;
    
        end
        if casenum == 1
            % at 8-20-07 11:50,
            r2ans = [5897.954130507     5791.046114526     6682.733686585];
            v2ans = [  -4.393910234     4.576816355     1.482423676];
    
            dut1  =  -0.1639883;
            dat   = 33;
            xp    =  0.210428 * convrt;
            yp    =  0.286899 * convrt;
            lod   =  0.0;
            timezone= 0;
            terms = 0;
            ddpsi = 0.0;
            ddeps = 0.0;
            obs1 = 6;
            obs2 = 10;
            obs3 = 14;
        end
        if casenum == 2
            % at 8-20-12 11:48:28.000 center time,
            %             year  =  2012;
            %             mon   =   8;
            %             day   =  20;
            %             hr    =  11;
            %             min   =  55;
            %             sec   =  28.0000;
            dut1  =  -0.6096413;
            dat   = 35; % leap second in July 2012
            xp    =  0.137495 * convrt;
            yp    =  0.342416 * convrt;
            lod   =  0.0;
            timezone= 0;
            terms = 0;
            ddpsi = 0.0;
            ddeps = 0.0;
        end
        % set eop parameters for new cases
        if casenum == 3
            dut1  =  0; %-0.6096413;
            dat   = 34; % leap second in July 2012
            xp    =  0.0; %137495 * convrt;
            yp    =  0.0; %342416 * convrt;
            lod   =  0.0;
            timezone= 0;
            terms = 0;
            ddpsi = 0.0;
            ddeps = 0.0;
    
            obs1 = 1; %62
            obs2 = 2; %70  72
            obs3 = 3; %74 102
        end
    
        if casenum == 4
            dut1  =  -0.1069721;
            dat   = 37;
            xp    =  0.148470 * convrt;
            yp    =  0.246564 * convrt;
            lod   =  0.0;
            timezone= 0;
            terms = 0;
            ddpsi = 0.0; %-0.113370; units
            ddeps = 0.0; %-0.007262;
    
            obs1 = 1; %62
            obs2 = 2; %70  72
            obs3 = 3; %74 102
            %2021 11 12 59530  0.148470  0.246564 -0.1069721  0.0001709 -0.113370 -0.007262  0.000222 -0.000052  37
        end
    
        switch(casenum)
            case 0
                answ="# 0  ans a = 12246.023  e = 0.2000  i = 40.00  330.000          km";
            case 1
                answ="# 1  old Book      2.0  ---- COEs: a= 12756.274 e= 0.20000 i= 62.3000   40.00   20.00 Topocentric case";
            case 2
                answ="# 2  ESC pg 282   1.5 ---- COEs:  a= 8306.247 e= 0.16419 i= 32.8780  136.53  203.95 Topocentric case -2519.36 8068.22 -2664.48 -5.9319 0.07621 2.60215";
                diffsites = 'y';
            otherwise
                fprintf('Invalid option\n' );
        end
    
    
        % ------ read all the data in and process
        %numobs = 3; % just get # of rows
    
    
        %             utc = 0.0;
        %             ut1 = utc+dut1;
        %             tai = utc+dat;
        %             tt  = tai+32.184;
        %             [jdut1, jdut1frac] = jday(year,mon,day,hr,min,ut1);
        %             [jdtt, jdttfrac]  = jday(year,mon,day,hr,min,tt);
        %             ttt   =  (jdtt-2451545.0)/36525.0;
        %             fprintf(1,'year %5i ',year);
        %             fprintf(1,'mon %4i ',mon);
        %             fprintf(1,'day %3i ',day);
        %             fprintf(1,'hr %3i:%2i:%8.6f\n',hr,min,second );
        %             fprintf(1,'dut1 %8.6f s',dut1);
        %             fprintf(1,' dat %3i s',dat);
        %             fprintf(1,' xp %8.6f "',xp);
        %             fprintf(1,' yp %8.6f "',yp);
        %             fprintf(1,' lod %8.6f s\n',lod);
    
    
        %load data into x y z arrays
        yeararr = filedat(:,3);
        monarr  = filedat(:,2);
        dayarr  = filedat(:,1);
        hrarr   = filedat(:,4);
        minarr  = filedat(:,5);
        secarr  = filedat(:,6);
        latarr  = filedat(:,7)/rad;
        lonarr  = filedat(:,8)/rad;
        altarr  = filedat(:,9); % km
        rtascarr  = filedat(:,10)/rad; % rad
        declarr   = filedat(:,11)/rad; % rad
    
        for j = 1:numobs
            [obsrecarr(j,1).jd,obsrecarr(j,1).jdf] = jday(yeararr(j),monarr(j),dayarr(j),hrarr(j),minarr(j),secarr(j));
            obsrecarr(j,1).latgd = latarr(j);  % assumes the same sensor site
            obsrecarr(j,1).lon = lonarr(j);
            obsrecarr(j,1).alt = altarr(j);
    
            utc = secarr(j);
            ut1 = utc + dut1;
            tai = utc + dat;
            tt  = tai + 32.184;
            [jdut1, jdut1frac] = jday(yeararr(j), monarr(j), dayarr(j), hrarr(j), minarr(j),ut1);
            [jdtt, jdttfrac]  = jday(yeararr(j), monarr(j), dayarr(j), hrarr(j), minarr(j),tt);
            ttt   =  (jdtt-2451545.0)/36525.0;
    
            [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac ] ...
                = convtime ( yeararr(j), monarr(j), dayarr(j), hrarr(j), minarr(j), secarr(j), 0, dut1, dat );
    
            fprintf(1,'time %i %i %i %i %i %i %f ttt %f %12.9f %f %12.9f %f %12.9f %f \n',j, ...
                yeararr(j), monarr(j), dayarr(j), hrarr(j), minarr(j), ...
                ut1, tt, ttt, jdut1, jdut1frac, jdtt, jdttfrac );
    
            [obsrecarr(j,1).rs,obsrecarr(j,1).vs] = site ( latarr(j),lonarr(j),altarr(j) );
            obsrecarr(j,1).ttt = ttt;
            obsrecarr(j,1).jdut1 = jdut1;
            obsrecarr(j,1).jdut1frac = jdut1frac;
            obsrecarr(j,1).xp = xp;  % rad
            obsrecarr(j,1).yp = yp;
            obsrecarr(j,1).rtasc = rtascarr(j);
            obsrecarr(j,1).decl = declarr(j);
        end
    
        rtasc1 = obsrecarr(obs1,1).rtasc;
        rtasc2 = obsrecarr(obs2,1).rtasc;
        rtasc3 = obsrecarr(obs3,1).rtasc;
        decl1 = obsrecarr(obs1,1).decl;
        decl2 = obsrecarr(obs2,1).decl;
        decl3 = obsrecarr(obs3,1).decl;
        jd1 = obsrecarr(obs1,1).jd;
        jdf1 = obsrecarr(obs1,1).jdf;
        jd2 = obsrecarr(obs2,1).jd;
        jdf2 = obsrecarr(obs2,1).jdf;
        jd3 = obsrecarr(obs3,1).jd;
        jdf3 = obsrecarr(obs3,1).jdf;
        rs1 = obsrecarr(obs1,1).rs;  % ecef
        vs1 = obsrecarr(obs1,1).vs;
        rs2 = obsrecarr(obs2,1).rs;
        vs2 = obsrecarr(obs2,1).vs;
        rs3 = obsrecarr(obs3,1).rs;
        vs3 = obsrecarr(obs3,1).vs;
    
        [year,mon,day,hr,min,second] = invjday(obsrecarr(obs1,1).jd, obsrecarr(obs1,1).jdf);
    
    
        % -------------- convert each site vector from ecef to eci -----------------
        a = [0;0;0];   % dummy acceleration variable for the ecef2eci routine
        [year,mon,day,hr,min,sec] = invjday(jd1,jdf1);
        [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac ] ...
            = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        [rsite1,vseci,aeci] = ecef2eci(rs1,vs1,a,iau80arr, ttt, jdut1+jdut1frac,lod,xp,yp,2,ddpsi,ddeps);
    
        [year,mon,day,hr,min,sec] = invjday(jd2,jdf2);
        [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac ] ...
            = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        [rsite2,vseci,aeci] = ecef2eci(rs2,vs2,a,iau80arr, ttt, jdut1+jdut1frac,lod,xp,yp,2,ddpsi,ddeps);
    
        [year,mon,day,hr,min,sec] = invjday(jd3,jdf3);
        [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac ] ...
            = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        [rsite3,vseci,aeci] = ecef2eci(rs3,vs3,a, iau80arr, ttt, jdut1+jdut1frac, lod, xp, yp, 2, ddpsi, ddeps );

    
        % ---------------------- run the angles-only routine ------------------
        if typerun == 'l' || typerun == 'a'
            [r2,v2] = anglesl( decl1,decl2,decl3,rtasc1,rtasc2, ...
                rtasc3,jd1,jdf1,jd2,jdf2,jd3,jdf3, diffsites, rsite1,rsite2,rsite3 );
            %                          rtasc3,jd1,jd2,jd3, rs1,rs2,rs3, re, mu, tu );
            processtype = 'anglesl';
            % -------------- write out answer --------------
            fprintf(1,'\ninputs: \n\n');
            [latgc,latgd,lon,alt] = ecef2ll ( rs1 ); % need to use ecef one!!
            fprintf(1,'Site obs1 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n', rsite1, latgd*rad, lon*rad, alt*1000 );
            [latgc,latgd,lon,alt] = ecef2ll ( rs2 );
            fprintf(1,'Site obs2 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n', rsite2, latgd*rad, lon*rad, alt*1000 );
            [latgc,latgd,lon,alt] = ecef2ll ( rs3 );
            fprintf(1,'Site obs3 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n', rsite3, latgd*rad, lon*rad, alt*1000 );
            [year,mon,day,hr,min,sec] = invjday ( jd1, jdf1 );
            fprintf(1,'obs#1 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n', year,mon,day,hr,min,sec, rtasc1*rad, decl1*rad );
            [year,mon,day,hr,min,sec] = invjday ( jd2, jdf2 );
            fprintf(1,'obs#2 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n', year,mon,day,hr,min,sec, rtasc2*rad, decl2*rad );
            [year,mon,day,hr,min,sec] = invjday ( jd3, jdf3 );
            fprintf(1,'Obs#3 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n', year,mon,day,hr,min,sec, rtasc3*rad, decl3*rad );
    
            fprintf(1,'\nsolution by %s \n\n', processtype);
            fprintf(1,'r2     %11.7f   %11.7f  %11.7f er    %11.7f  %11.7f  %11.7f km \n',r2/re, r2);
            %        fprintf(1,'r2 ans %11.7f   %11.7f  %11.7f er    %11.7f  %11.7f  %11.7f km \n',r2ans/re, r2ans);
    
            fprintf(1,'v2     %11.7f   %11.7f  %11.7f er/tu %11.7f  %11.7f  %11.7f km/s\n',v2/velkmps, v2);
            %        fprintf(1,'v2 ans %11.7f   %11.7f  %11.7f er/tu %11.7f  %11.7f  %11.7f km/s\n',v2ans/velkmps, v2ans);
    
            [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r2,v2, re, mu);
            fprintf(1,'         p km          a km         ecc       incl deg     raan deg    argp deg     nu deg      m deg  \n');
            fprintf(1,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f \n',...
                p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    
            fprintf(1,'%s \n', answ );
        end  % laplace
    
        if typerun == 'g' || typerun == 'a'
            [r2,v2] = anglesg( decl1,decl2,decl3,rtasc1,rtasc2, ...
                rtasc3,jd1,jdf1,jd2,jdf2,jd3,jdf3, rsite1,rsite2,rsite3 );
            processtype = 'anglesg';
            % -------------- write out answer --------------
            fprintf(1,'\ninputs: \n\n');
            [latgc,latgd,lon,alt] = ecef2ll ( rs1 ); % need to use ecef one!!
            fprintf(1,'Site obs1 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n', rsite1, latgd*rad, lon*rad, alt*1000 );
            [latgc,latgd,lon,alt] = ecef2ll ( rs2 );
            fprintf(1,'Site obs2 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n', rsite2, latgd*rad, lon*rad, alt*1000 );
            [latgc,latgd,lon,alt] = ecef2ll ( rs3 );
            fprintf(1,'Site obs3 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n', rsite3, latgd*rad, lon*rad, alt*1000 );
            [year,mon,day,hr,min,sec] = invjday ( jd1, jdf1 );
            fprintf(1,'obs#1 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n', year,mon,day,hr,min,sec, rtasc1*rad, decl1*rad );
            [year,mon,day,hr,min,sec] = invjday ( jd2, jdf2 );
            fprintf(1,'obs#2 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n', year,mon,day,hr,min,sec, rtasc2*rad, decl2*rad );
            [year,mon,day,hr,min,sec] = invjday ( jd3, jdf3 );
            fprintf(1,'Obs#3 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n', year,mon,day,hr,min,sec, rtasc3*rad, decl3*rad );
    
            fprintf(1,'\nsolution by %s \n\n', processtype);
            fprintf(1,'r2     %11.7f   %11.7f  %11.7f er    %11.7f  %11.7f  %11.7f km \n',r2/re, r2);
            %        fprintf(1,'r2 ans %11.7f   %11.7f  %11.7f er    %11.7f  %11.7f  %11.7f km \n',r2ans/re, r2ans);
    
            fprintf(1,'v2     %11.7f   %11.7f  %11.7f er/tu %11.7f  %11.7f  %11.7f km/s\n',v2/velkmps, v2);
            %        fprintf(1,'v2 ans %11.7f   %11.7f  %11.7f er/tu %11.7f  %11.7f  %11.7f km/s\n',v2ans/velkmps, v2ans);
    
            [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r2,v2, re, mu);
            fprintf(1,'         p km          a km         ecc       incl deg     raan deg    argp deg     nu deg      m deg  \n');
            fprintf(1,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f \n',...
                p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    
            fprintf(1,'%s \n', answ );
        end  % gauss
        %pause;
    
        if typerun == 'd' || typerun == 'a'
            pctchg = 0.05;
    
            bigr2x = getGaussRoot(decl1,decl2,decl3,rtasc1,rtasc2, ...
                rtasc3,jd1,jdf1,jd2,jdf2,jd3,jdf3, rsite1,rsite2,rsite3);
            if (bigr2x < 0.0 || bigr2x > 50000.0)
                bigr2x = 35000.0;  % simply set this to about GEO, allowing for less than that too.
                % use Laplace guess
            end
            bigr2x
            rng1 = bigr2x;
            rng2 = bigr2x * 1.02;
    
            [r2,v2] = anglesdr( decl1,decl2,decl3,rtasc1,rtasc2, ...
                rtasc3,jd1,jdf1,jd2,jdf2,jd3,jdf3, diffsites,rsite1,rsite2,rsite3, rng1, rng2, pctchg);
            % -------------- write out answer --------------
            fprintf(1,'\ninputs: \n\n');
            [latgc,latgd,lon,alt] = ecef2ll ( rs2 ); % need to use ecef one!!
            fprintf(1,'rSite eci %11.7f %11.7f %11.7f km  \n', rsite1 );
            fprintf(1,'Site obs1 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n', rsite1, latgd*rad, lon*rad, alt*1000 );
            [latgc,latgd,lon,alt] = ecef2ll ( rs2 );
            fprintf(1,'Site obs2 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n', rsite2, latgd*rad, lon*rad, alt*1000 );
            [latgc,latgd,lon,alt] = ecef2ll ( rs3 );
            fprintf(1,'Site obs3 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n', rsite3, latgd*rad, lon*rad, alt*1000 );
            [year,mon,day,hr,min,sec] = invjday ( jd1, jdf1 );
            fprintf(1,'obs#1 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n', year,mon,day,hr,min,sec, rtasc1*rad, decl1*rad );
            [year,mon,day,hr,min,sec] = invjday ( jd2, jdf2 );
            fprintf(1,'obs#2 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n', year,mon,day,hr,min,sec, rtasc2*rad, decl2*rad );
            [year,mon,day,hr,min,sec] = invjday ( jd3, jdf3 );
            fprintf(1,'Obs#3 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n', year,mon,day,hr,min,sec, rtasc3*rad, decl3*rad );
    
            fprintf(1,'r2     %11.7f   %11.7f  %11.7f er    %11.7f  %11.7f  %11.7f km \n',r2/re, r2);
            %        fprintf(1,'r2 ans %11.7f   %11.7f  %11.7f er    %11.7f  %11.7f  %11.7f km \n',r2ans/re, r2ans);
    
            fprintf(1,'v2     %11.7f   %11.7f  %11.7f er/tu %11.7f  %11.7f  %11.7f km/s\n',v2/velkmps, v2);
            %        fprintf(1,'v2 ans %11.7f   %11.7f  %11.7f er/tu %11.7f  %11.7f  %11.7f km/s\n',v2ans/velkmps, v2ans);
    
            [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r2,v2);
            fprintf(1,'         p km          a km         ecc       incl deg     raan deg    argp deg     nu deg      m deg  \n');
            fprintf(1,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f \n',...
                p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    
            fprintf(1,'%s \n', answ );
    
            fprintf(1,'=============================================== \n');
            pause;
        end  % double -r
    
    end  % through all the cases
    
    
