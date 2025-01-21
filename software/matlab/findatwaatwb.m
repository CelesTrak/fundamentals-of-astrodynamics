    % -----------------------------------------------------------------------------
    %
    %                           procedure findatwbatwa
    %
    % this procedure finds the a and b matrices for the differential correction
    %   problem.  remember that it isn't critical for the propagations to use
    %   the highest fidelity techniques because we're only trying to find the
    %   "slope". k is an index that allows us to do multiple rows at once. it's
    %   used for both the b and a matrix calculations.
    %
    %  algorithm     : find the a and b matrices by accumulation to reduce matrix sizes
    %                  calculate the matrix combinations
    %                  atw is found without matrix operations to avoid large matrices
    %
    %  author        : david vallado                  719-573-2600    6 aug 2008
    %
    %  inputs          description                                range / units
    %    firstob     - number of observations
    %    lastob      - number of observations
    %    statesize   - size of state                                6 , 7
    %    percentchg  - amount to modify the vectors
    %                  by in finite differencing
    %    deltaamtchg - tolerance for small value in
    %                  finite differencing                         0.0000001
    %    whichconst  - parameter for sgp4 constants                wgs72, wgs721, wgs84
    %    satrec      - structure of satellite parameters for TLE
    %    obsrecfile    - array of records containing:
    %                  senum, jd, rsvec, obstype,
    %                  rng, az, el, drng, daz, del,
    %                  trtasc, tdecl data
    %    statetype   - type of elements (equinoctial, etc)         'e', 't'
    %    scalef      - scale factor to limit the size of the state vector
    %    xnom        - state vector                                varied
    %
    %  outputs       :
    %    atwa        - atwa matrix
    %    atwb        - atwb matrix
    %    atw         - atw matrix
    %    b           - b matrix, the residuals
    %    drng2       - range residual squared
    %    daz2        - azimuth residual squared
    %    del2        - elevation residual squared
    %    ddrng2      - range rate residual squared
    %    ddaz2       - azimuth rate residual squared
    %    ddel2       - elevation rate residual squared
    %    dtrtasc2    - topocentric right ascension residual squared
    %    dtdecl2     - topocentric declination residual squared
    %    dx2         - x position residual squared
    %    dy2         - y position residual squared
    %    dz2         - z position residual squared
    %    dxdot2      - xdot position residual squared
    %    dydot2      - ydot position residual squared
    %    dzdot2      - zdot position residual squared
    %
    %  locals        :
    %    rnom        - nom position vector at epoch                km
    %    vnom        - nom velocity vector at epoch                km/s
    %    a           - a matrix
    %    indobs      -
    %    at          -
    %    w1          -
    %    w2          -
    %    w3          -
    %    lst         -
    %    gst         -
    %    dtsec        -
    %    deltaamt    -
    %    rngpert     -
    %    azpert      - modified azimuth                            -2pi to 2pi
    %    elpert      - modified azimuth                            -pi/2 to pi/2
    %    drng        -
    %    daz         -
    %    del         -
    %    error       -
    %
    %  coupling      :
    %    findsenptr  - find sensor data
    %    rv_razel    - find r and v given range, az, el, and rates
    %    rv_tradec   - find r and v given topocentric rtasc and decl
    %
    %  references    :
    %    vallado       2013, 753-765
    % --------------------------------------------------------------------------- */

    function [atwa, atwb, atw, b, drng2, daz2, del2 ] = findatwaatwb(iau80arr, firstobs, lastobs, obsrecarr,  statesize, percentchg, deltaamtchg, xnom);

    % --------------------- initialize parameters ------------------
    if obsrecarr(1).obstype == 0
        indobs = 1;
    elseif obsrecarr(1).obstype == 2
        indobs = 3;
    else
        indobs = 2;
    end

    drng2 = 0.0;
    daz2 = 0.0;
    del2 = 0.0;
    dtrtasc2 = 0.0;
    dtdecl2 = 0.0;
    atwaacc = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
    atwbacc = [0; 0; 0; 0; 0; 0;];

    % ------------- reset these since they will accumulate ---------
    % zero out matrices
    for r = 1: statesize
        for c = 1: statesize
            atwa(r,c) = 0.0;
        end
        atwb(r,1) = 0.0;
    end

    % ------------------- loop through all the observations ------------------
    for obsktr = firstobs: lastobs
        currobsrec = obsrecarr(obsktr);
        %printf( "ob %2i rsecef0 %11.5f jd %11.5f dtmin %8.3f hr %3i rng %11.7f \n",
        %         i, currobsrec.rsecef[0], currobsrec.jd, currobsrec.dtmin, currobsrec.hr, currobsrec.rng );  % ritrf

        % ----------------- determine sensor characteristics -----------------
        %         rs(1) = currobsrec.rsecef(1);
        %         rs(2) = currobsrec.rsecef(2);
        %         rs(3) = currobsrec.rsecef(3);

        % temporary sensor for now
        %  getsensorparams( currobsrec.sennum, currsenrec );


        % --------- propagate the nominal vector to the epoch time -----------
        %         sgp4 (whichconst, satrec,  currobsrec.dtmin, rteme, vteme);
        dtsec = (obsrecarr(obsktr).time - obsrecarr(1).time + obsrecarr(obsktr).timef - obsrecarr(1).timef) * 86400.0;  % s
        rnom(1) = xnom(1,1);
        rnom(2) = xnom(2,1);
        rnom(3) = xnom(3,1);
        vnom(1) = xnom(4,1);
        vnom(2) = xnom(5,1);
        vnom(3) = xnom(6,1);
        [reci1, veci1] =  kepler ( rnom, vnom, dtsec );
        %[reci1, veci1] =  pkepler ( rnom, vnom, dtsec, 0.0, 0.0 );
        aeci = [0;0;0];
        [recef1, vecef1, aecef1] = eci2ecef(reci1', veci1', aeci, iau80arr, currobsrec.ttt,currobsrec.jdut1,0.0,currobsrec.xp,currobsrec.yp,2,0.0,0.0 );
        [rsecef,vsecef] = site ( currobsrec.latgd,currobsrec.lon,currobsrec.alt);

        % ------------------------- find b matrix ----------------------------
        if currobsrec.obstype ~= 3
            [rngnom,aznom,elnom,drngnom,daznom,delnom] = rv2razel ( recef1,vecef1, currobsrec.latgd,currobsrec.lon,currobsrec.alt);
        else
            [rngnom,trtascnom,tdeclnom,drngnom,dtrtascnom,dtdeclnom] = rv2tradec(recef1, vecef1, rsecef, vsecef);
        end

        switch (currobsrec.obstype)
            case 0
                b(1,1) = currobsrec.rng - rngnom;
            case 1
                b(1,1) = currobsrec.az  - aznom;
                %fix for 0-360...
                if (abs(b(1,1)) > pi)
                    b(1,1) = b(1,1) - sgn(b(1,1)) * 2.0 * pi;
                end
                b(2,1) = currobsrec.el  - elnom;
            case 2
                b(1,1) = currobsrec.rng  - rngnom;
                b(2,1) = currobsrec.az   - aznom;
                % fix for 0-360...
                if abs(b(2,1)) > pi
                    b(2,1) = b(2,1) - sign(b(2,1)) * 2.0 * pi;
                end
                b(3,1) = currobsrec.el  - elnom;
            case 3
                b(1,1) = currobsrec.trtasc - trtascnom;
                % fix for 0-360...
                if (abs(b(1,1)) > pi)
                    b(1,1) = b(1,1) - sign(b(1,1)) * 2.0 * pi;
                end
                b(2,1) = currobsrec.tdecl  - tdeclnom;
        end  
        %printf( "rnom %11.5f %11.5f %11.5f %8.3f %8.3f %8.3f %8.3f \n",
        %              rteme[0], rteme[1], rteme[2], rngnom, aznom*rad, elnom*rad, currobsrec.rng );  % ritrf


        % ------------------------ find a matrix -----------------------------
        % ------------- reset the perturbed vector to the nominal ------------
        %  satrecp = satrec;
        xnomp = xnom;

        % ----- perturb each element in the state (elements or vectors) ------
        for j= 1: statesize
            [deltaamt, xnomp] = finitediff(j, percentchg, deltaamtchg, xnom);

            %   sgp4 (whichconst, satrecp,  currobsrec.dtmin, r3, v3);
            dtsec = (obsrecarr(obsktr).time - obsrecarr(1).time + obsrecarr(obsktr).timef - obsrecarr(1).timef) * 86400;
            rnomp(1) = xnomp(1,1);
            rnomp(2) = xnomp(2,1);
            rnomp(3) = xnomp(3,1);
            vnomp(1) = xnomp(4,1);
            vnomp(2) = xnomp(5,1);
            vnomp(3) = xnomp(6,1);
            [reci3, veci3] =  kepler ( rnomp, vnomp, dtsec );
            [recef3, vecef3, aece31] = eci2ecef(reci3', veci3', aeci, iau80arr, currobsrec.ttt,currobsrec.jdut1,0.0,currobsrec.xp,currobsrec.yp,2,0.0,0.0 );
            %[reci3, veci3] =  pkepler ( rnomp, vnomp, dtsec, 0.0, 0.0 );

            % teme to itrf if observation type
            %  if (currobsrec.obstype ~= 4)
            %      mfme = currobsrec.hr(j) * 60 + currobsrec.min(j) + currobsrec.sec(j)/60.0;
            %      findeopparam ( currobsrec.jd, mfme, interp, eoparr, jdeopstart,
            %                     dut1, dat, lod, xp, yp, ddpsi, ddeps,
            %                     iaudx, dy, icrsx, y, s, deltapsi, deltaeps );
            %      [ut1, tut1, jdut1, utc, tai, tt, ttt, jdtt, tdb, ttdb, jdtdb ] ...
            %               = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
            %      iau76fk5_itrf_teme( ritrf, vitrf, aitrf, eFrom, r3, v3, ateme, ttt, xp, yp, jdut1, lod, trans );
            %  end % if obstype

            if (currobsrec.obstype == 3)
                [trrpert,trtascpert,tdeclpert,tdrrpert,tdrtascpert,tddeclpert] = rv2tradec(recef3, vecef3, rsecef, vsecef);
            else
                if (currobsrec.obstype == 4)
                    bstarpert = satrec.bstar*(1.0 + percentchg);
                else  % currobsrec.obstype = 0 or 1 or 2;
                    [rngpert,azpert,elpert,drhopert,dazpert,delpert] = rv2razel ( recef3,vecef3, currobsrec.latgd,currobsrec.lon,currobsrec.alt);
                    %        fprintf(1,'rnom  %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f km \n',reci1, rngnom,aznom,elnom, deltaamt );
                    %        fprintf(1,'rpert %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f km \n',reci3,rngpert,azpert,elpert, dtsec );
                end
            end
            switch currobsrec.obstype
                case 0
                    a(1,j) = (rngpert - rngnom) / deltaamt;
                case 1
                    a(1,j) = (azpert  - aznom)  / deltaamt;
                    a(2,j) = (elpert  - elnom)  / deltaamt;
                case 2
                    a(1,j) = (rngpert - rngnom) / deltaamt;
                    a(2,j) = (azpert  - aznom)  / deltaamt;
                    a(3,j) = (elpert  - elnom)  / deltaamt;
                case 3
                    a(1,j) = (trtascpert - trtascnom) / deltaamt;
                    a(2,j) = (tdeclpert  - tdeclnom)  / deltaamt;
            end % case

            %printf( "rpert %11.5f %11.5f %11.5f %8.3f %8.3f %8.3f \n",
            %              r3[0], r3[1], r3[2], rngpert, azpert*rad, elpert*rad );  % ritrf


            % ----------------- reset the modified vector --------------------
            %satrecp = satrec;
            xnomp = xnom;
        end  % for j = 0 to statesize

        % ----------------- now form the matrix combinations -----------------
        at = a';

        % ------------------------- assign weights ---------------------------
        switch (currobsrec.obstype)
            case 0
                w1 = 1.0 / (currsenrec.noiserng * currsenrec.noiserng);
                rng2 = rng2 + b(1,1) * b(1,1) * w1;
            case 1
                w1 = 1.0 / (currsenrec.noiseaz * currsenrec.noiseaz);
                w2 = 1.0 / (currsenrec.noiseel * currsenrec.noiseel);
                daz2 = daz2 + b(1,1) * b(1,1) * w1;
                del2 = del2 + b(2,1) * b(2,1) * w2;
            case 2
                w1 = 1.0 / (currobsrec.noiserng * currobsrec.noiserng);
                w2 = 1.0 / (currobsrec.noiseaz * currobsrec.noiseaz);
                w3 = 1.0 / (currobsrec.noiseel * currobsrec.noiseel);
                drng2 = drng2 + b(1,1) * b(1,1) * w1;
                daz2  = daz2  + b(2,1) * b(2,1) * w2;
                del2  = del2  + b(3,1) * b(3,1) * w3;
            case 3
                w1 = 1.0 / (currobsrec.noisetrtasc * currobsrec.noisetrtasc);
                w2 = 1.0 / (currobsrec.noisetdecl * currobsrec.noisetdecl);
                dtrtasc2 = dtrtasc2 + b(1,1) * b(1,1) * w1;
                dtdecl2  = dtdecl2  + b(2,1) * b(2,1) * w2;
        end % case

        for rowc = 1: statesize
            for colc = 1: indobs
                weight = 1.0;
                switch (colc)
                    case 1
                        weight= w1;
                    case 2
                        weight= w2;
                    case 3
                        weight= w3;
                    case 4
                        weight= w4;
                    case 5
                        weight= w5;
                    case 6
                        weight= w6;
                    case 7
                        weight= w7;
                end  % case
                atw(rowc,colc) = at(rowc,colc) * weight;
            end  % for colc
        end % for rowc

        % ----------------- find the atwa / atwb matrices --------------------
        atwaacc = atw * a; % don't add atwaacc twice!!
        %matmult( atw,a, atwaacc, statesize,indobs,statesize );
        atwbacc = atw * b;
        %matmult( atw,b, atwbacc, statesize,indobs,1 );

        % ------------------- accumulate the matricies -----------------------
        for r = 1: statesize
            for c = 1: statesize
                atwa(r,c) = atwaacc(r,c) + atwa(r,c);
            end
        end

        c = 1;
        for r = 1: statesize
            atwb(r,c) = atwbacc(r,c) + atwb(r,c);
        end % for obsktr through the observations

        %writeexpmat("atwa ",atwa,statesize,statesize);
        %writeexpmat("atwb ",atwb,statesize,1);
        %writeexpmat("a ",a,indobs,statesize);
        %writemat("b ",b,indobs,1);
    end  % for obsktr through the observations


