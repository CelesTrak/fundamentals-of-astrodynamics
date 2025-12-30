    % script testsgp4.m
    %
    % This script tests the SGP4 propagator.

    % Author:
    %   Jeff Beck
    %   beckja@alumni.lehigh.edu

    % Version Info:
    %   1.0 (051019) - Initial version from Vallado C++ version.
    %   1.0 (aug 14, 2006) - update for paper
    %   2.0 (apr 2, 2007) - update for manual operations
    %   3.0 (3 jul, 2008) - update for opsmode operation afspc or improved
    %   3.1 (2 dec, 2008) - fix tsince/1440.0 in jd update

    % sgp4fix consolidate call to getgravconst in sgp4init
    % these are set in sgp4init
    % global tumin mu radiusearthkm xke j2 j3 j4 j3oj2

    directory = 'd:\codes\library\matlab\';
    fprintf(1,'output directory set to: %s, change in testsgp4.m if needed \n', directory);

    % global opsmode

    % // ------------------------  implementation   --------------------------

    %   add operation smode for afspc (a) or improved (i)
    %opsmode= input('input opsmode afspc (a), improved i ','s');
    % this is the standard method of operation
    opsmode = 'a'; 
    
    %         //typerun = 'c' compare 1 year of full satcat data
    %         //typerun = 'v' verification run, requires modified elm file with
    %         //typerun = 'm' maunual operation- either mfe, epoch, or dayof yr
    %         //              start stop and delta times
    typerun = input('input type of run c, v, m: ','s');
    if (typerun == 'm')
        typeinput = input('input mfe, epoch (YMDHMS), or dayofyr approach, m,e,d: ','s');
    else
        typeinput = 'e';
    end

    % whichconst = input('input constants 721, (72), 84 ');
    % this is the standard method of operation
    whichconst = 72;
   
    rad = 180.0 / pi;

    % placeholder
    ateme = [0.0 0.0 0.0];

    %         // ---------------- setup files for operation ------------------
    %         // input 2-line element set file
    infilename = input('input elset filename: ','s');
    infile = fopen(infilename, 'r');
    if (infile == -1)
        fprintf(1,'Failed to open file: %s\n', infilename);
        return;
    end

    if (typerun == 'c')
        outfile = fopen(strcat(directory,'testsgp4all.out'), 'wt');
    else
        if (typerun == 'v')
            outfile = fopen(strcat(directory,'testsgp4.out'), 'wt');
        else
            outfile = fopen(strcat(directory,'tsgp4.out'), 'wt');
        end
    end
    if (outfile == -1)
        fprintf(1,'Failed to open outfile in %s\n', directory);
        return;
    end

    global idebug dbgfile
    mu         = 398600.8;      % km3/s2, WGS072 value

    
    
    % alternatively, if you want to find the 2 lines for a specific satellite from a
    % file, use the following:
    %sscnumstr = '1 43234U';
    %text = fileread(infilename);
    %TextAsCells = regexp(text, '\n', 'split');
    %mask = contains(TextAsCells, sscnumstr);
    %longstr1 = TextAsCells(mask);
    %sscnumstr = '2 43234U';
    %mask = contains(TextAsCells, sscnumstr);
    %longstr2 = TextAsCells(mask);
    % Note that Mat;ab requires some funny thinsg to get a string array...
    %[startmfe, stopmfe, deltamin, satrec] = twoline2rv( ...
    %   convertStringsToChars(string(longstr1)), convertStringsToChars(string(longstr2)), typerun, typeinput, opsmode, whichconst);


    
    
    % ----------------- test simple propagation -------------------
    while (~feof(infile))
        longstr1 = fgets(infile, 130);
        while ( (longstr1(1) == '#') && (feof(infile) == 0) )
            longstr1 = fgets(infile, 130);
        end

        if (feof(infile) == 0)

            longstr2 = fgets(infile, 130);


            % sgp4fix additional parameters to store from the TLE
            satrec.classification = 'U';
            satrec.intldesg = '        ';
            satrec.ephtype = 0;
            satrec.elnum   = 0;
            satrec.revnum  = 0;

            if idebug
                catno = strtrim(longstr1(3:7));
                dbgfile = fopen(strcat('sgp4test.dbg.',catno), 'wt');
                fprintf(dbgfile,'this is the debug output\n\n' );
            end
            %  convert the char string to sgp4 elements
            %  includes initialization of sgp4
            [startmfe, stopmfe, deltamin, satrec] = twoline2rv( ...
                longstr1, longstr2, typerun, typeinput, opsmode, whichconst);

            % sgp4fix note that the ephtype must be 0 for SGP4. SGP4-XP uses 4.
            if (satrec.ephtype == 0)

                fprintf(outfile, '%d xx\n', satrec.satnumStr);
                fprintf(1,' %d\n', satrec.satnumStr);

                % call the propagator to get the initial state vector value
                [satrec, rteme ,vteme] = sgp4 (satrec,  0.0);

                fprintf(outfile, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n',...
                    satrec.t,rteme(1),rteme(2),rteme(3),vteme(1),vteme(2),vteme(3));
                %            fprintf(outfile, ' %16.7f %16.7f %16.7f %16.7f %12.8f %12.8f %12.8f\n',...
                %                 satrec.t,ro(1),ro(2),ro(3),vo(1),vo(2),vo(3));
                %            fprintf(1, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n',...
                %                 satrec.t,ro(1),ro(2),ro(3),vo(1),vo(2),vo(3));

                tsince = startmfe;

                %                // check so the first value isn't written twice
                if ( abs(tsince) > 1.0e-8 )
                    tsince = tsince - deltamin;
                end

                % set some EOP parameter that stay constant
                conv = pi / (180.0*3600.0);
                eqeterms = 2;  % terms for equation of the equinoxes

                % loop to perform the propagation
                while ((tsince < stopmfe) && (satrec.error == 0))

                    tsince = tsince + deltamin;

                    if(tsince > stopmfe)
                        tsince = stopmfe;
                    end

                    [satrec, rteme, vteme] = sgp4 (satrec,  tsince);
                    if (satrec.error > 0)
                        fprintf(1,'# *** error: t:= %f *** code = %3i\n', tsince, satrec.error);
                    end

                    if (satrec.error == 0)
                        if ((typerun ~= 'v') && (typerun ~= 'c'))
                            jdutc = satrec.jdsatepoch;
                            jdutcfrac = satrec.jdsatepochf + tsince/1440.0;
                            if jdutcfrac < 0.0
                                jdutc = jdutc - 1.0;
                                jdutcfrac = jdutcfrac + 1.0;
                            end
                            [year,mon,day,hr,minute,sec] = invjday ( jdutc, jdutcfrac );

                            fprintf(outfile,...
                                ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f %5i%3i%3i %2i:%2i:%9.6f  \n',...
                                tsince,rteme(1),rteme(2),rteme(3),vteme(1),vteme(2),vteme(3),year,mon,day,hr,minute,sec );
                            % %16.8f%16.8f%16.8%12.9f%12.9f%12.9f and \r

                            % demo getting lat lon
                            % set the EOP parameters. normally this would be at each time instant in the loop, but
                            % we'll keep it constant here for illustration
                            lod = 0.0015563; % sec
                            xp   = -0.140682 * conv;  % polar motion values in rad from arcsec
                            yp   =  0.333309 * conv;
                            dut1 = -0.4399619;  % sec
                            jdut1 = jdutc;
                            jdut1frac = jdutcfrac + dut1/86400.0;
                            ttt = (jdut1-2451545.0)/36525.0;
                            [recef, vecef, aecef] = teme2ecef( rteme', vteme', ateme', ttt, jdut1, lod, xp, yp, eqeterms );
                            %[latgc,latgd,lon,hellp] = ecef2ll ( recef );
                            %fprintf(1,' lat lon %11.7f  %11.7f \n', latgc*rad, lon*rad);
                        else
                            jdutc = satrec.jdsatepoch;
                            jdutcfrac = satrec.jdsatepochf + tsince/1440.0;
                            if jdutcfrac < 0.0
                                jdutc = jdutc - 1.0;
                                jdutcfrac = jdutcfrac + 1.0;
                            end
                            [year,mon,day,hr,minute,sec] = invjday ( jdutc, jdutcfrac );

                            fprintf(outfile, ' %16.8f %16.8f %16.8f %16.8f %12.8f %12.8f %12.8f',...
                                tsince,rteme(1),rteme(2),rteme(3),vteme(1),vteme(2),vteme(3));
                            %                        fprintf(1, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f \n',...
                            %                            tsince,ro(1),ro(2),ro(3),vo(1),vo(2),vo(3));

                            [p,a,ecc,incl,node,argp,nu,m,arglat,truelon,lonper ] = rv2coe (rteme,vteme);

                            fprintf(outfile, ' %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f %5i%3i%3i %2i:%2i:%9.6f \n',...
                                a, ecc, incl*rad, node*rad, argp*rad, nu*rad, m*rad,year,mon,day,hr,minute,sec );
                            % demo getting lat lon.
                            % set some EOP parameters. normally this would be at each time instant in the loop, but
                            % we'll keep it constant here for illustration
                            lod = 0.0015563; % sec
                            xp   = -0.140682 * conv;  % polar motion values in rad from arcsec
                            yp   =  0.333309 * conv;
                            dut1 = -0.4399619;  % sec
                            jdut1 = jdutc;
                            jdut1frac = jdutcfrac + dut1/86400.0;
                            ttt = (jdut1-2451545.0)/36525.0;
                            [recef, vecef, aecef] = teme2ecef( rteme', vteme', ateme', ttt, jdut1, lod, xp, yp, eqeterms );
                            %[latgc,latgd,lon,hellp] = ecef2ll ( recef );
                            %fprintf(1,' lat lon %11.7f  %11.7f \n', latgc*rad, lon*rad);
                        end
                    end % if satrec.error == 0

                end % while propagating the orbit

                if (idebug && (dbgfile ~= -1))
                    fclose(dbgfile);
                end
            end


        end % if not eof

    end % while through the input file

    fclose(infile);
    fclose(outfile);

    clear all;  % clear out the variables to test direct approach
    % sgp4fix demonstrate method of running SGP4 directly from orbital element values
    %1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813
    %2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656
    deg2rad  =   pi / 180.0;         %   0.0174532925199433
    xpdotp   =  1440.0 / (2.0 *pi);  % 229.1831180523293

    whichconst = 72;
    opsmode = 'a';
    satrec.satnumStr = '08195';
    satrec.jdsatepoch = 2453911.0;
    satrec.jdsatepochf = 0.8321544402;
    satrec.no_kozai = 2.00491383;
    satrec.ecco = 0.6877146;
    satrec.inclo = 64.1586;
    satrec.nodeo = 279.0717;
    satrec.argpo = 264.7651 ;
    satrec.mo = 20.2257;
    satrec.nddot = 0.00000e0;
    satrec.bstar = 0.11873e-3;
    satrec.ndot = 0.00000099;
    satrec.elnum = 813;
    satrec.revnum = 22565;
    satrec.classification = 'U';
    satrec.intldesg = '75081A';
    satrec.ephtype = 0;

    % convert units and initialize
    satrec.no_kozai = satrec.no_kozai / xpdotp; %* rad/min
    satrec.ndot  = satrec.ndot  / (xpdotp*1440.0);  %* ? * minperday
    satrec.nddot = satrec.nddot / (xpdotp*1440.0*1440);
    satrec.inclo = satrec.inclo  * deg2rad;
    satrec.nodeo = satrec.nodeo  * deg2rad;
    satrec.argpo = satrec.argpo  * deg2rad;
    satrec.mo    = satrec.mo     * deg2rad;

    % set start/stop times for propagation
    startmfe =     0.0;
    stopmfe  =  2880.0;
    deltamin =   120.0;

    [satrec] = sgp4init( whichconst, opsmode, satrec, satrec.jdsatepoch+satrec.jdsatepochf-2433281.5, satrec.bstar, ...
        satrec.ndot, satrec.nddot, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no_kozai, ...
        satrec.nodeo);

    tsince = startmfe;
    while ((tsince < stopmfe) && (satrec.error == 0))
        tsince = tsince + deltamin;

        if (tsince > stopmfe)
            tsince = stopmfe;
        end

        [satrec, rteme, vteme] = sgp4 (satrec,  tsince);

        jdutc = satrec.jdsatepoch;
        jdutcfrac = satrec.jdsatepochf + tsince/1440.0;
        if jdutcfrac < 0.0
            jdutc = jdutc - 1.0;
            jdutcfrac = jdutcfrac + 1.0;
        end
        [year,mon,day,hr,minute,sec] = invjday ( jdutc,jdutcfrac );

        fprintf( 1,' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f \n',...
            tsince,rteme(1),rteme(2),rteme(3),vteme(1),vteme(2),vteme(3));
    end % while propagating the orbit

