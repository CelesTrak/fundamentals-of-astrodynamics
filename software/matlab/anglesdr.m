%  ------------------------------------------------------------------------------
%
%                           procedure anglesdoubler
%
%  this procedure solves the problem of orbit determination using three
%    optical sightings. the solution procedure uses the double-r technique.
%    the important thing is the input of the initial guesses of the range which
%    may be easiest from the solution of the gauss 8th order poly.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    tdecl1       - declination #1                               rad
%    tdecl2       - declination #2                               rad
%    tdecl3       - declination #3                               rad
%    trtasc1      - right ascension #1                           rad
%    trtasc2      - right ascension #2                           rad
%    trtasc3      - right ascension #3                           rad
%    trtasc3      - right ascension #3                           rad
%    jd1, jdf1    - julian date of 1st sighting                  days from 4713 bc
%    jd2, jdf2    - julian date of 2nd sighting                  days from 4713 bc
%    jd3, jdf3    - julian date of 3rd sighting                  days from 4713 bc
%    rs1          - eci site position vector #1                  km
%    rs2          - eci site position vector #2                  km
%    rs3          - eci site position vector #3                  km
%    magr1in      - initial estimate
%    magr2in      - initial estimate
%
%  outputs         :
%    r2           -  position vector at t2                       km
%    v2           -  velocity vector at t2                       km / s
%    errstr       - output results for debugging
%
%  locals         :
%    los1         - line of sight vector for 1st
%    los2         - line of sight vector for 2nd
%    los3         - line of sight vector for 3rd
%    tau          - taylor expansion series about
%                   tau ( t - to )
%    tausqr       - tau squared
%    t21t23       - (t2-t1) * (t2-t3)
%    t31t32       - (t3-t1) * (t3-t2)
%    rho          - range from site to sat at t2                 km
%    rhodot       -
%    earthrate    - velocity of earth rotation
%
%  coupling       :
%    mag          - magnitude of a vector
%    matmult      - multiply two matrices together
%    angle        - angle between two vectors
%
%  references     :
%    vallado       2022, 450, alg 53, ex 7-2
%
% [r2, v2] = anglesdr(decl1, decl2, decl3, rtasc1, rtasc2, ...
%        rtasc3, jd1, jdf1, jd2, jdf2, jd3, jdf3, diffsites, rseci1, rseci2, rseci3, rng1, rng2, pctchg);
% ------------------------------------------------------------------------------

function [r2, v2] = anglesdr(decl1, decl2, decl3, rtasc1, rtasc2, ...
        rtasc3, jd1, jdf1, jd2, jdf2, jd3, jdf3, diffsites, rseci1, rseci2, rseci3, rng1, rng2, pctchg)

    % -------------------------  implementation   -------------------------
    constastro;

    %show = 'y';
    show = 'n';

    % for sun
    %re = 149597870.0;
    %mu = 1.32712428e11;
    rad = 180.0 / pi;

    magr1in = rng1;  % 2.0*re;
    magr2in = rng2;  % 2.04*re;
    %direct  = 'y';

    tol    = 1e-8*re;   % km
    tol = 0.1; % km
    %pctchg = 0.005;

    % subtract dates and convert fraction of day to seconds
    tau12 = (jd1 - jd2) * 86400.0 + (jdf1 - jdf2) * 86400.0; % days to sec
    tau13 = (jd3 - jd1) * 86400.0 + (jdf3 - jdf1) * 86400.0;
    tau32 = (jd3 - jd2) * 86400.0 + (jdf3 - jdf2) * 86400.0;

    period = 86400.0;
    n12 = floor( abs(tau12 / period) );
    n13 = floor( abs(tau13 / period) );
    n23 = floor( abs((tau12 + tau32) / period) );

    % form line of sight vectors
    los1 = [cos(decl1)*cos(rtasc1) cos(decl1)*sin(rtasc1) sin(decl1)]';
    los2 = [cos(decl2)*cos(rtasc2) cos(decl2)*sin(rtasc2) sin(decl2)]';
    los3 = [cos(decl3)*cos(rtasc3) cos(decl3)*sin(rtasc3) sin(decl3)]';

    % --------- now we're ready to start the actual double r algorithm ---------
    magr1old  = 99999.9;
    magr2old  = 99999.9;
    magrsite1 = mag(rseci1);
    magrsite2 = mag(rseci2);
    magrsite3 = mag(rseci3);

    % take away negatives because escobal defines rs opposite
    cc1 = 2.0*dot(los1,rseci1);
    cc2 = 2.0*dot(los2,rseci2);
    ktr = 0;

    oldqr = 100000000.0;
    newqr = 9000000.0;
    % main loop to get three values of the double-r for processing
    while (abs(magr1in-magr1old) > tol | abs(magr2in-magr2old) > tol && ktr <= 15 ...
            && newqr < oldqr)

        oldqr = newqr;
        ktr = ktr + 1;
        if show == 'y'
            fprintf(1,'%2i ',ktr);
        end
        magr1o = magr1in;
        magr2o = magr2in;

        [r2, r3, f1, f2, q1, magr1, magr2, a, deltae32] = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,...
            los1, los2, los3, rseci1, rseci2, rseci3, tau12, tau32, n12, n13, n23);

        % check intermediate status
        f  = 1.0 - a/magr2*(1.0-cos(deltae32));
        g  = tau32 - sqrt(a^3/mu)*(deltae32-sin(deltae32));
        v2 = (r3 - f*r2)/g;
        [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r2,v2);
        if show == 'y'
            fprintf(1,'coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
                p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
        end

        % -------------- re-calculate f1 and f2 with r1 = r1 + delta r1
        deltar1 = pctchg * magr1o;
        magr1in = magr1o + deltar1;
        magr2in = magr2o;
        [r2, r3, f1delr1, f2delr1, q2, magr1, magr2, a, deltae32] = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,...
            los1, los2, los3, rseci1, rseci2, rseci3, tau12, tau32, n12, n13, n23);

        pf1pr1 = (f1delr1-f1)/deltar1;
        pf2pr1 = (f2delr1-f2)/deltar1;

        % ----------------  re-calculate f1 and f2 with r2 = r2 + delta r2
        magr1in = magr1o;
        deltar2 = pctchg * magr2o;
        magr2in = magr2o + deltar2;
        [r2, r3, f1delr2, f2delr2, q3, magr1, magr2, a, deltae32] = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,...
            los1, los2, los3, rseci1, rseci2, rseci3, tau12, tau32, n12, n13, n23);

        pf1pr2 = (f1delr2-f1)/deltar2;
        pf2pr2 = (f2delr2-f2)/deltar2;

        %       f  = 1.0 - a/magr2*(1.0-cos(deltae32));
        %       g  = t3 - sqrt(a^3/mu)*(deltae32-sin(deltae32));
        %       v2 = (r3 - f*r2)/g;
        %       [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r2,v2);
        %       fprintf(1,'coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        %              p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );


        % ------------ now calculate an update
        magr2in = magr2o;

        delta  = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
        delta1 = pf2pr2*f1 - pf1pr2*f2;
        delta2 = pf1pr1*f2 - pf2pr1*f1;

        if (abs(delta) < 0.0000001)
            deltar1 = -delta1;
            deltar2 = -delta2;
        else
            deltar1 = -delta1/delta;
            deltar2 = -delta2/delta;
        end

        magr1old = magr1in;
        magr2old = magr2in;

        %  may need to limit the amount of the correction
        chkamt = 0.15;
        if abs(deltar1 / magr1in) > chkamt
            deltar1 = magr1in * chkamt * sign(deltar1);
            if show == 'y'
                fprintf(1,'%11.7f \n',deltar1 );
            end
            %         deltar1 = sign(deltar1)*magr1in*pctchg;
        end
        if abs(deltar2 / magr2in) > chkamt
            deltar2 = magr2in * chkamt * sign(deltar2);
            if show == 'y'
                fprintf(1,'%11.7f \n',deltar2 );
            end
            %         deltar2 = sign(deltar2)*magr2in*pctchg;
        end

        magr1old = magr1in;
        magr2old = magr2in;

        magr1in = magr1in + deltar1;
        magr2in = magr2in + deltar2;

        if show == 'y'
            fprintf(1,'qs %11.7f  %11.7f  %11.7f \n',q1,q2,q3);
            fprintf(1,'magr1o %11.7f delr1 %11.7f magr1 %11.7f %11.7f  \n',magr1o,deltar1,magr1in,magr1old);
            fprintf(1,'magr2o %11.7f delr2 %11.7f magr2 %11.7f %11.7f  \n',magr2o,deltar2,magr2in,magr2old);
        end

        newqr = sqrt(q1*q1+q2*q2+q3*q3);

        %       f  = 1.0 - a/magr2*(1.0-cos(deltae32));
        %       g  = t3 - sqrt(a^3/mu)*(deltae32-sin(deltae32));
        %       v2 = (r3 - f*r2)/g;
        %       [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (r2,v2);
        %       fprintf(1,'coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        %              p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );

        if show == 'y'
            fprintf(1,'=============================================== \n');
        end

        pctchg = pctchg * 0.5;

        %pause;
    end

    % needed to get the r2 set properly since the last one was moving r2
    [r2, r3, f1, f2, q1, magr1, magr2, a, deltae32] = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in,...
        los1, los2, los3, rseci1, rseci2, rseci3, tau12, tau32, n12, n13, n23);

    f  = 1.0 - a/magr2*(1.0-cos(deltae32));
    g  = tau32 - sqrt(a^3/mu)*(deltae32-sin(deltae32));
    v2 = (r3 - f*r2)/g;

end