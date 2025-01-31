    %     -----------------------------------------------------------------
    %
    %                              Ex4_1.m
    %
    %  this file demonstrates example 4-1.
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

    fprintf(1,'--------------- book angle conversion tests ----------------------------\n' );
    latgd = 39.007/rad;
    lon = -104.883/rad;
    alt = 2.19456; % km

    year = 2004;
    mon  =   5;
    day  =  14;
    hr   =    12;
    min  =   0;
    sec  =   0.00;

    year = 1994;
    mon  =   5;
    day  =  14;
    hr   =    13;
    min  =   11;
    sec  =   20.59856;

    dut1 =  0.0;
    dat  = 32;
    xp   =  0.0;
    yp   =  0.0;
    lod  =  0.0;
    timezone = 0;
    order =  106;
    terms = 2;

    fileLoc = 'D:\Codes\LIBRARY\DataLib\';
    [iau80arr] = iau80in(fileLoc);

    % , tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac
    [ ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, ...
      tdb, ttdb, jdtdb, jdtdbfrac ] ...
      = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
    fprintf(1,'ut1 %8.6f tut1 %16.12f jdut1 %18.11f\n',ut1,tut1,jdut1 );
    fprintf(1,'utc %8.6f\n',utc );
    fprintf(1,'tai %8.6f\n',tai );
    fprintf(1,'tt  %8.6f ttt  %16.12f jdtt  %18.11f\n',tt,ttt,jdtt );
    fprintf(1,'tdb %8.6f ttdb %16.12f jdtdb %18.11f\n',tdb,ttdb,jdtdb );

    [lst, gst] = lstime(lon, jdut1);
    fprintf(1,'lst %11.7f gst %11.7f \n',lst*rad, gst*rad );

    for i = 1:2
        if i == 1
            fprintf(1,'\n-------- Neptune test baseline test \n' );
            reci = [1752246215.0; -3759563433.0; -1577568105.0];
            veci = [-18.324; 18.332; 7.777 ];
            aeci = [0.001;0.002;0.003];
            reci = reci;
            veci = veci;
            rr    =  29.664361*149597870.0;
            rtasc = 294.98914583/rad;
            decl  = -20.8234944/rad;
            % old book value                drr   = (149598023.0*(29.649616 - 29.664361))/86400.0
            drr   = (149597870.0*(29.649616 - 29.664361))/86400.0
            drtasc= -0.00000012244/rad;
            ddecl = -0.00000001794/rad;
            [reci,veci] = radec2rv(rr,rtasc,decl,drr,drtasc,ddecl)
        end
        if i == 2
            fprintf(1,'\n-------- closer test baseline test \n' );
            rr    =  12756.0;
            rtasc = 294.98914583/rad;
            decl  = -20.8234944/rad;
            drr   = 6.798514;
            drtasc= -0.00000012244/rad;
            ddecl = -0.00000001794/rad;
            [reci,veci] = radec2rv(rr,rtasc,decl,drr,drtasc,ddecl);
        end
        % geoc
        fprintf(1,'r    %14.7f%14.7f%14.7f',reci );
        fprintf(1,' v %14.9f%14.9f%14.9f\n',veci );

        [recef, vecef, aecef] = eci2ecef(reci, veci, aeci, iau80arr, ttt, jdut1, lod, xp, yp, 2, ddpsi, ddeps );
        fprintf(1,'ITRF rev      IAU-76/FK5  %15.11f  %15.11f  %15.11f %15.11f  %15.11f  %15.11f\n', ...
            recef(1), recef(2), recef(3), vecef(1), vecef(2), vecef(3));


        [rr,rtasc,decl,drr,drtasc,ddecl] = rv2radec( reci,veci );
        fprintf(1,'            rho km       rtasc deg     decl deg      drho km/s      drtasc deg/s   ddecl deg/s \n' );
        if rtasc < 0.0
            rtasc = rtasc + twopi;
        end
        fprintf(1,'radec  %14.7f %14.7f %14.7f',rr,rtasc*rad,decl*rad );
        fprintf(1,' %14.7f %14.12f %14.12f\n',drr,drtasc*rad,ddecl*rad );

        [reci,veci] = radec2rv(rr,rtasc,decl,drr,drtasc,ddecl);
        fprintf(1,'r    %14.7f %14.7f %14.7f',reci );
        fprintf(1,' v %14.9f %14.9f %14.9f\n',veci );

        % topoc
        % ----------------- get site vector in ecef -------------------
        [rsecef, vsecef] = site ( latgd, lon, alt );
        % -------------------- convert ecef to eci --------------------
        a = [0;0;0];
        [rseci, vseci, aeci] = ecef2eci(rsecef, vsecef, a, iau80arr, ttt, jdut1, lod, xp, yp, 2, ddpsi, ddeps);

        [trr, trtasc, tdecl, ddrr, tdrtasc, tddecl] = rv2tradec ( recef, vecef, rsecef, vsecef );
        fprintf(1,'           trho km      trtasc deg    tdecl deg     tdrho km/s     tdrtasc deg/s  tddecl deg/s \n' );
        if trtasc < 0.0
            trtasc = trtasc + twopi;
        end
        fprintf(1,'tradec  %14.7f %14.7f %14.7f',trr,trtasc*rad,tdecl*rad );
        fprintf(1,' %14.7f %14.12f %14.12f \n',tdrr,tdrtasc*rad,tddecl*rad );

        [r1, v1] = tradec2rv (trr, trtasc, tdecl, tdrr, tdrtasc, tddecl, rseci', vseci');
        fprintf(1,'tradec r    %14.7f%14.7f%14.7f',r1 );
        fprintf(1,' v %14.9f%14.9f%14.9f \n',v1 );

        %horiz
        [rho, az, el, drho, daz, del] = rv2razel(recef, vecef, latgd, lon, alt );
        if az < 0.0
            az = az + twopi;
        end
        fprintf(1,'rvraz   %14.7f %14.7f %14.7f',rho,az*rad,el*rad );
        fprintf(1,' %14.7f %14.12f %14.12f\n',drho,daz*rad,del*rad );

        [recef, vecef] = razel2rv(latgd, lon, alt, rho, az, el, drho, daz, del);
        [rseci, vseci, aeci] = ecef2eci(rsecef, vsecef, a, iau80arr, ttt, jdut1, lod, xp, yp, 2, ddpsi, ddeps);
        fprintf(1,'r    %14.7f %14.7f %14.7f',reci );
        fprintf(1,' v %14.9f %14.9f %14.9f\n',veci );


        % ecl lat lon
        [rr,elon,elat,drr,delon,delat] = rv2elatlon( reci,veci );
        fprintf(1,'            rho km        elon deg     elat deg      drho km/s       delon deg/s   delat deg/s \n' );
        if elon < 0.0
            elon = elon + twopi;
        end
        fprintf(1,'ell      %14.7f %14.7f %14.7f',rr,elon*rad,elat*rad );
        fprintf(1,' %14.7f %14.12f %14.12f\n',drr,delon*rad,delat*rad );

        [reci,veci] = elatlon2rv(rr,elon,elat,drr,delon,delat);
        fprintf(1,'r    %14.7f %14.7f %14.7f',reci );
        fprintf(1,' v %14.9f %14.9f %14.9f\n',veci );
    end % for
    

    
    

    
    
