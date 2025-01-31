% ------------------------------------------------------------------------------
%
%                           function eq2rv
%
%  this function finds the classical orbital elements given the equinoctial
%    elements.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    a           - semimajor axis                             km
%    af          - component of ecc vector
%    ag          - component of ecc vector
%    chi         - component of node vector in eqw
%    psi         - component of node vector in eqw
%    meanlon     - mean longitude                             rad
%    fr          - retrograde factor, neg if incl > 179 deg 1, -1
%
%  outputs       :
%    r           - position vector                            km
%    v           - velocity vector                            km/s
%
%  locals        :
%    n           - mean motion                                rad
%    temp        - temporary variable
%    p           - semilatus rectum                           km
%    ecc         - eccentricity
%    incl        - inclination                                0.0  to pi rad
%    raan        - rtasc of ascending node                    0.0  to 2pi rad
%    argp        - argument of perigee                        0.0  to 2pi rad
%    nu          - true anomaly                               0.0  to 2pi rad
%    m           - mean anomaly                               0.0  to 2pi rad
%    arglat      - argument of latitude      (ci)             0.0  to 2pi rad
%    truelon     - true longitude            (ce)             0.0  to 2pi rad
%    lonper      - longitude of periapsis    (ee)             0.0  to 2pi rad
%
%  coupling      :
%
%  references    :
%    vallado 2022 : 110
%
% [r, v] = eq2rv( a, af, ag, chi, psi, meanlonM, fr);
% ------------------------------------------------------------------------------

function [r, v] = eq2rv( a, af, ag, chi, psi, meanlonM, fr)
    % -------------------------  implementation   -----------------
    constmath;
    constastro;

    arglat  = 999999.1;
    lonper  = 999999.1;
    truelon = 999999.1;

    coe = true;  % pick coe or vector approaches

    if (coe == true)

        % ---- if n is input ----
        %a = (mu/n^2)^(1.0/3.0);

        ecc = sqrt (af^2 + ag^2);
        p = a * (1.0 - ecc*ecc);
        incl = pi*((1.0 - fr)*0.5) + 2.0*fr*atan( sqrt(chi^2 + psi^2) );
        omega = atan2( chi, psi);
        argp = atan2( ag,af ) - fr*omega;

        if ( ecc < small )
            % ----------------  circular equatorial  ------------------
            if (incl<small) || ( abs(incl-pi)< small )
                %argp = 0.0;
                %omega= 0.0;
                truelon = omega;
            else
                % --------------  circular inclined  ------------------
                %argp= 0.0;
                arglat = argp;
            end
        else
            % ---------------  elliptical equatorial  -----------------
            if ( ( incl<small) || (abs(incl-pi)<small) )
                %argp = lonper;
                lonper = omega;  % ok
            end
        end

        m = meanlonM - fr*omega - argp;
        m = rem (m + twopi,twopi);

        [e0, nu] = newtonm ( ecc, m );

        % ----------  fix for elliptical equatorial orbits ------------
        if ( ecc < small )
            % ----------------  circular equatorial  ------------------
            if (incl<small) || ( abs(incl-pi)< small )
                argp    = undefined;
                omega   = undefined;
                truelon = nu;
            else
                % --------------  circular inclined  ------------------
                argp  = undefined;
                arglat= nu - fr*omega;
            end
            nu   = undefined;
        else
            % ---------------  elliptical equatorial  -----------------
            if ( ( incl < small) || (abs(incl-pi) < small) )
                lonper = argp;           % ok
                argp    = undefined;
                omega   = undefined;
            end
        end

        % -------- now convert back to position and velocity vectors
        [r, v] = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper);

    else
        p0 = 1.0 / (1.0 + chi^2 + psi^2);
        fe = p0 * (1.0 - chi^2 + psi^2);  % 2nd one is minus???? no, + seems correct
        fq = p0 * 2.0 * chi * psi;
        fw = p0 * -2.0 * fr*chi;
        fvec = [fe, fq, fw];
        ge = p0 * 2.0 * fr*chi * psi;
        gq = p0 * fr*(1.0 + chi^2 - psi^2);
        gw = p0 * 2.0 * psi;
        gvec = [ge, gq, gw];
        we = p0 * 2.0 * chi;
        wq = p0 * -2.0 * psi;
        ww = p0 * fr*(1.0 - chi^2 - psi^2);
        wvec = [we, wq, ww];  % same


        F0 = meanlonM;
        numiter = 25;
        ktr= 1;
        F1 = F0 - (F0 + ag*cos(F0) - af*sin(F0) - meanlonM) / (1.0 - ag*sin(F0) - af*cos(F0));
        while (( abs(F1-F0) > small ) && ( ktr <= numiter ))
            ktr = ktr + 1;
            F0= F1;
            F1 = F0 - (F0 + ag*cos(F0) - af*sin(F0) - meanlonM)/(1.0 - ag*sin(F0) - af*cos(F0));
            %            fprintf(1,'iters %7i  %11.7f  %11.7f \n', ktr, F0, F1);
        end;
        F = F1;
        F = rem(F+twopi,twopi);

        n  = sqrt(mu/(a*a*a));

        temp = (af^2 + ag^2);
        if abs(temp) > 1.0
            temp = 1.0;
        end
        b = 1.0 / (1.0 + sqrt(1.0 - temp));

        sinL = ((1.0 - af^2*b)*sin(F) + ag*af*b*cos(F) - ag) / (1.0 - ag*sin(F) - af*cos(F));
        cosL = ((1.0 - ag^2*b)*cos(F) + ag*af*b*sin(F) - af) / (1.0 - ag*sin(F) - af*cos(F));
        L = atan2(sinL, cosL);
        meanlonNu = L;
        meanlonNu = rem(meanlonNu + twopi, twopi);

        rr = a*(1.0 - ag^2 - af^2) / (1.0 + ag*sin(L) + af*cos(L));
        rr1 = a*(1.0 - ag*sin(F) - af * cos(F));

        % coordinates in equinoctial space
        X = a*((1.0 - ag^2 * b) * cos(F) + af*ag*b*sin(F) - af);
        Y = a*((1.0 - af^2 * b) * sin(F) + af*ag*b*cos(F) - ag);
        XD = -n*a*(ag + sinL) / (sqrt(1.0 - af^2 - ag^2));
        YD =  n*a*(af + cosL) / (sqrt(1.0 - af^2 - ag^2));
        % alt forms all are the same now
        %XD = n*a^2/magr * ( af*ag*b*cos(F) - (1.0 - ag^2*b)*sin(F) );
        %YD = n*a^2/magr * ( (1.0 - af^2*b)*cos(F) - af*ag*b*sin(F) );


        %        fprintf(1,'eq2rv fe %11.7f %11.7f %11.7f ge %11.7f  %11.7f %11.7f \n X %11.7f Y %11.7f b %11.7f sF %11.7f cF %11.7f \n', fe, fq, fw, ge, gq, gw, X, Y, b, sin(F), cos(F));

        %         r = a*(1.0-af^2-ag^2) / (1.0 + ag*sinL + af*cosL);
        %         r = a*(1.0 - ag*sin(F) - af*cos(F));
        %         r = magr
        %         X = r*cosL;
        %         Y = r*sinL;

        r = X*fvec + Y*gvec;
        v = XD*fvec + YD*gvec;
        %   fprintf(1,'eq2rv F %11.7f  L %11.7f \n', F*rad, L*rad);
    end


    % test for eqw axes
    %     r = [6524.834000000,  6862.875000000,  6448.296000000];
    %     v = [4.9013270000,    5.5337560000,   -1.9763410000];
    %
    %     incl = 87.8691262/rad;
    %     raan = 227.89826    /rad;
    %
    %
    %     [outvec] = rot3( r, raan );
    %     [outvec1] = rot1( outvec, incl );
    %     [ans] = rot3( outvec1, -raan )
    %    ans =  1.0e+04 * 1.113447759019948   0.269748810750719  0.000000005640130 correct

end