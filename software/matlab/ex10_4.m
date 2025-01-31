    %     -----------------------------------------------------------------
    %
    %                              Ex10_4.m
    %
    %  this file demonstrates example 10-4.
    %
    %                          companion code for
    %             fundamentals of astrodyanmics and applications
    %                                 2022
    %                            by david vallado
    %
    %     (h)               email davallado@gmail.com
    %
    %     *****************************************************************
    %
    %  current :
    %            19 feb 19  david vallado
    %                         update for new constants
    %  changes :
    %            16 oct 23  david vallado
    %                         recheck oepration
    %
    %     *****************************************************************

    constmath;
    constastro;

    fprintf(1,'problem 10-4 -------------------------\n');

    % read in GEOS data
    %     fid = fopen('d:\Codes\LIBRARY\Matlab\geos6.inp');
    %     obsarr = textscan(fid,'%d %d %d %d %d %d %d %d %f %f %f %f ', 'Delimiter', ' ', 'Headerlines', 1); % these are read in as cells
    %     fclose(fid);
    %
    %     yeararr = 1 * obsarr{4end(:);
    %     monarr  = 1 * obsarr{5end(:);
    %     dayarr  = 1 * obsarr{6end(:);
    %     hrarr   = 1 * obsarr{7end(:);
    %     minarr  = 1 * obsarr{8end(:);
    %     secarr  = 1 * obsarr{9end(:);
    %     rngarr  = obsarr{10end(:)/1000.0; % m
    %     azarr   = obsarr{11end(:)/rad; % deg
    %     elarr   = obsarr{12end(:)/rad; % deg


    filedat =load('d:\Codes\LIBRARY\Matlab\geos6a.inp');
    numobs = size(filedat,1); % just get # of rows
  %  numobs = 10;

    %load data into x y z arrays
    yeararr = filedat(:,4);
    monarr  = filedat(:,5);
    dayarr  = filedat(:,6);
    hrarr   = filedat(:,7);
    minarr  = filedat(:,8);
    secarr  = filedat(:,9);
    rngarr  = filedat(:,10); % km
    azarr   = filedat(:,11)/rad; % deg
    elarr  = filedat(:,12)/rad; % deg

    convrt = pi / (3600.0*180.0);
    latgd = 21.572056/rad;   % AAS paper values
    lon   = -158.266578/rad;
    alt   = 0.3002;  % km
    dat   = 29;  % assume the EOP values stay the same throughout the interval
    dut1  = 0.3261068;  % s
    lod   = 0.0;
    xp    = -0.11554 * convrt; % arcsec to rad
    yp    =  0.48187 * convrt;
    [rs,vs] = site ( latgd,lon,alt );

    for j = 1:numobs  
        [obsrecarr(j,1).time, obsrecarr(j,1).timef] = jday(yeararr(j),monarr(j),dayarr(j),hrarr(j),minarr(j),secarr(j));
        obsrecarr(j,1).latgd = latgd;  % assumes the same sensor site
        obsrecarr(j,1).lon = lon;
        obsrecarr(j,1).alt = alt;
        
        [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] ...
            = convtime ( yeararr(j), monarr(j), dayarr(j), hrarr(j), minarr(j), secarr(j), 0, dut1, dat );
        obsrecarr(j,1).ttt = ttt;
        obsrecarr(j,1).jdut1 = jdut1 + jdut1frac;  % just keep together for now
        obsrecarr(j,1).xp = xp;  % rad
        obsrecarr(j,1).yp = yp;
        obsrecarr(j,1).noiserng = 0.0925;  % km
        obsrecarr(j,1).noiseaz = 0.0224/rad;   % rad
        obsrecarr(j,1).noiseel = 0.0139/rad;   % rad
        obsrecarr(j,1).obstype = 2;  % range az el
        obsrecarr(j,1).rsecef = rs;
        obsrecarr(j,1).rng = rngarr(j) - 0.08;  % km subtract known biases off these values
        obsrecarr(j,1).az = azarr(j) - 0.0081/rad;
        obsrecarr(j,1).el = elarr(j) - 0.0045/rad;
    end

    firstobs = 1;
    lastobs = 10;

    % -------  form nominal vector in eci
    % simply average each vector moved back to the nominal time
    reci = [0; 0; 0];
    veci = [0; 0; 0];
    for obsktr = firstobs+1: lastobs-1
        currobsrec = obsrecarr(obsktr-1);
        [re1, ve1] = razel2rv(currobsrec.latgd,currobsrec.lon,currobsrec.alt, ...
            currobsrec.rng,currobsrec.az,currobsrec.el, 0.0, 0.0, 0.0);

        currobsrec = obsrecarr(obsktr);
        [re2, ve2] = razel2rv(currobsrec.latgd,currobsrec.lon,currobsrec.alt, ...
            currobsrec.rng,currobsrec.az,currobsrec.el, 0.0, 0.0, 0.0);

        currobsrec = obsrecarr(obsktr+1);
        [re3, ve3] = razel2rv(currobsrec.latgd,currobsrec.lon,currobsrec.alt, ...
            currobsrec.rng,currobsrec.az,currobsrec.el, 0.0, 0.0, 0.0);
        % far apart vectors
        % [ve2, theta,theta1,copa, error] = gibbs( re1,re2,re3);
        % closely spaced vectors
        [ve2, theta,theta1,copa, error] = hgibbs( re1,re2,re3,obsrecarr(obsktr-1).jdut1,obsrecarr(obsktr).jdut1,obsrecarr(obsktr+1).jdut1 );

        % move back to 1st time, not central time
        dtsec = (obsrecarr(obsktr).time + obsrecarr(obsktr).timef - obsrecarr(1).time - obsrecarr(1).timef) * 86400.0;  % s
        [reci1, veci1] =  kepler ( re2, ve2, -dtsec );
        %        [reci1, veci1] =  pkepler ( re2, ve2, -dtsec, 0.0, 0.0 );
        fprintf(1,'re %16.8f %16.8f %16.8f km ',reci1 );
        fprintf(1,'ve %16.8f %16.8f %16.8f %11.7f km \n',veci1,dtsec );

        reci = reci + reci1;
        veci = veci + veci1;
    end
    reci = reci/(lastobs-firstobs-1);
    veci = veci/(lastobs-firstobs-1);

    fprintf(1,'rnom %16.8f %16.8f %16.8f km ',reci );
    fprintf(1,'vnom %16.8f %16.8f %16.8f km \n',veci );

    numobs

    jdepoch = obsrecarr(1,1).time + obsrecarr(1,1).timef;
    % use a vector that's farther off to see iterations take effect
    reci = [5975.2904  2568.6400  3120.5845];
    veci = [3.983846  -2.071159  -5.917095];

    fprintf(1,'input: \n' );
    fprintf(1,'rnom %16.8f %16.8f %16.8f km ',reci );
    fprintf(1,'vnom %16.8f %16.8f %16.8f km \n',veci );

    xnom(1,1) = reci(1);
    xnom(2,1) = reci(2);
    xnom(3,1) = reci(3);
    xnom(4,1) = veci(1);
    xnom(5,1) = veci(2);
    xnom(6,1) = veci(3);

    % set parameters for finite differencing
    percentchg = 0.010; % 0.001
    deltaamtchg = 0.01; % 0.0000001
    for j = 1:5  % 5 iterations for now

        % ---- accumulate obs and assemble matrices
        [atwa, atwb, atw, b, drng2, daz2, del2] = findatwaatwb(firstobs, lastobs, obsrecarr, 6, percentchg, deltaamtchg, xnom);

        if j == 1
            fprintf(1,'atwa = \n' );
            fprintf(1,'%10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n',atwa );
            fprintf(1,'atwb = \n' );
            fprintf(1,'%10.1f \n%10.1f \n%10.1f \n%10.1f \n%10.1f \n%10.1f\n',atwb );
        end
        
        % matrix inverse
        deltax = inv(atwa)*atwb;
        fprintf(1,'dx %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f \n',deltax );

        xnom = xnom + deltax;
        r1(1) = xnom(1,1);
        r1(2) = xnom(2,1);
        r1(3) = xnom(3,1);
        v1(1) = xnom(4,1);
        v1(2) = xnom(5,1);
        v1(3) = xnom(6,1);

        % answer in km and km/s
        fprintf(1,'output: \n' );
        fprintf(1,'r1 %16.8f %16.8f %16.8f km ',r1 );
        fprintf(1,'v1 %16.8f %16.8f %16.8f km \n',v1 );

        switch (obsrecarr(1).obstype)
            %            0 : sigmanew= SQRT( DDRng2 / NumWork );
            %            1 : sigmanew= SQRT( (DAz2 + DEl2) / NumWork );
            case 2
                sigmanew = sqrt( (drng2 + daz2 + del2) / (lastobs-firstobs) );
                %            3 : sigmanew= SQRT( (DRng2 + DAz2 + DEl2 + DDRng2) / NumWork );
                %            4 : sigmanew= SQRT( (DRng2 + DAz2 + DEl2 + DDRng2 + DDAz2 + DDEl2)
                %                                 / NumWork );
                %            5 : sigmanew= SQRT( (DTRtAsc2 + DTDecl2) / NumWork );
                %            6 : sigmanew= SQRT( DRng2 / NumWork );
        end % case
        fprintf(1,'rms  %16.8f  \n', sigmanew );

    end % through iterations

    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r1,v1);
    fprintf(1,'          p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n');
    fprintf(1,'coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad, ...
        arglat*rad,truelon*rad,lonper*rad );


    rans = [5748.6070  2679.9324  3442.7902];  % from book and gtds
    vans = [4.330460  -1.922874  -5.726562];

    r2shans = [5746.469300   2679.944362   3442.416160 ]; % odtk short answer w bad intiial
    v2shans = [4.33993195   -1.91866466   -5.72186789];

    r3ans = [5748.398438   2680.135601   3442.296315]; % odtk answer from all data in LS w bad initial rnom
    v3ans = [4.32743417   -1.92014180   -5.72605888]; % odtk answer from all data in LS

    dr = r2shans - r1;
    fprintf(1,'diff odtk short bad init  %16.8f  %16.8f  %16.8f  %16.8f  km \n', dr, mag(dr) );

    dr = rans - r1;
    fprintf(1,'diff gtds w bad init      %16.8f  %16.8f  %16.8f  %16.8f  km \n', dr, mag(dr) );

    dr = r3ans - r1;
    fprintf(1,'diff odtk long w bad init %16.8f  %16.8f  %16.8f  %16.8f  km \n', dr, mag(dr) );

    p = inv(atwa);
    fprintf(1,'covariance p = \n' );
    fprintf(1,'%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e \n',p );

    fprintf(1,'cov r, 1 sigma  %16.8f  %16.8f  %16.8f  km \n', sqrt(p(1,1))*1000, sqrt(p(2,2))*1000, sqrt(p(3,3))*1000 );

    [eigenaxes,d] = eig(p);
    eigenvalues = sqrt(d)*1000;  % in m
    fprintf(1,'eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n', eigenaxes(1:6) );
    fprintf(1,'eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n', eigenaxes(7:12) );
    fprintf(1,'eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n', eigenaxes(13:18) );
    fprintf(1,'eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n', eigenaxes(19:24) );
    fprintf(1,'eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n', eigenaxes(25:30) );
    fprintf(1,'eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n', eigenaxes(31:36) );
    fprintf(1,'eigenvalues   %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  m \n', eigenvalues(1,1), eigenvalues(2,2), eigenvalues(3,3), ...
        eigenvalues(4,4), eigenvalues(5,5), eigenvalues(6,6) );
    fprintf(1,'mag eigenvalues  %16.8f m  \n', sqrt(eigenvalues(1,1)^2+ eigenvalues(2,2)^2+ eigenvalues(3,3)^2)  );
    fprintf(1,'mag eigenvalues  %16.8f m/s \n', sqrt(eigenvalues(4,4)^2+ eigenvalues(5,5)^2+ eigenvalues(6,6)^2)  );

    rms = sqrt(b*b'/4);
    rmsold = rms;
    %     fprintf(1,'0 dx %11.7f  %11.7f ans %11.7f  %11.7f rms %11.7f \n',ans(1), ans(2), alpha, beta, rms);

    
    
