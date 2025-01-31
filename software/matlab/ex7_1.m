    %     -----------------------------------------------------------------
    %
    %                              Ex7_1.m
    %
    %  this file demonstrates example 7-1.
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
    %             7 jun 07  david vallado
    %                         original
    %  changes :
    %            13 feb 07  david vallado
    %                         original baseline
    %
    %     *****************************************************************
    
    constmath;
    
    % ---- site
    fprintf(1,'\n-------- site test \n' );
    latgd = 39.007/rad;
    lon = -104.883/rad;
    alt = 2.1870;
    [jd,jdfrac]  = jday(1995,5,20, 3,17,2.0);
    [rs,vs] = site ( latgd,lon,alt );
    fprintf(1,'site %14.7f%14.7f%14.7f%14.7f%14.7f%14.7f\n',rs,vs );
    
    fprintf(1,'--------------- razel tests ----------------------------\n' );
    latgd =  39.007/rad;
    lon   = -104.883/rad;
    alt   =  2.187;
    rho   =  604.68;  % km
    az    = 205.6/rad;
    el    = 30.7/rad;
    drho  = 2.08;
    daz   = 0.15/rad;
    del   = 0.17/rad;
    year  = 1995;
    mon   =   5;
    day   =  20;
    hr    =   3;
    min   =  17;
    sec   =   2.0;
    dut1  =  0.0;
    dat   = 29;
    xp    =  0.0;
    yp    =  0.0;
    lod   =  0.0;
    timezone=0;
    terms = 0;
    ddpsi = 0.0;
    ddeps = 0.0;
    
    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    utc=sec;
    ut1=utc+dut1;
    tai=utc+dat;
    tt=tai+32.184;
    [jdut1,jdut1frac]= jday(year,mon,day,hr,min,ut1);
    [jdtt,jdttfrac]= jday(year,mon,day,hr,min,tt);
    ttt =  (jdtt-2451545.0)/36525.0;
    
    fprintf(1,'year %5i ',year);
    fprintf(1,'mon %4i ',mon);
    fprintf(1,'day %3i ',day);
    fprintf(1,'hr %3i:%2i:%8.6f\n',hr,min,sec );
    fprintf(1,'dut1 %8.6f s',dut1);
    fprintf(1,' dat %3i s',dat);
    fprintf(1,' xp %8.6f "',xp);
    fprintf(1,' yp %8.6f "',yp);
    fprintf(1,' lod %8.6f s\n',lod);
    fprintf(1,'           range km        az deg      el    deg     rngrt km/s      azrate deg/s  elrate deg/s\n');
    fprintf(1,'rvraz %14.7f%14.7f%14.7f%14.7f%14.7f%14.7f\n',rho,az*rad,el*rad,drho,daz*rad,del*rad );
    
    [recef, vecef] = razel2rv(latgd, lon, alt, rho, az, el, drho, daz, del);
    [reci, veci, aeci] = ecef2eci(rsecef, vsecef, a, iau80arr, ttt, jdut1, lod, xp, yp, 2, ddpsi, ddeps);
    
    fprintf(1,'r    %14.7f%14.7f%14.7f',reci );
    fprintf(1,' v %14.9f%14.9f%14.9f\n',veci );
    
    [rho, az, el, drho, daz, del] = rv2razel(recef, vecef, latgd, lon, alt );
    
    fprintf(1,'rvraz %14.7f%14.7f%14.7f%14.7f%14.7f%14.7f\n',rho,az*rad,el*rad,drho,daz*rad,del*rad );
    
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (reci, veci);
    fprintf(1,'          p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n');
    fprintf(1,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad, ...
        arglat*rad,truelon*rad,lonper*rad );
    

