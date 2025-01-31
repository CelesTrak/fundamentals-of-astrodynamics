% ----------------------------------------------------------------------------
%
%                           function rv2eq
%
% this function transforms a position and velocity vector into the flight
%    elements - latgc, lon, fpa, az, position and velocity magnitude.
%
%  author        : david vallado             davallado@gmail.com      20 jan 2025
%
%  inputs          description                              range / units
%    r           - eci position vector                        km
%    v           - eci velocity vector                        km/s
%
%  outputs       :
%    n           - mean motion                                rad
%    a           - semi major axis                            km
%    af          - component of ecc vector
%    ag          - component of ecc vector
%    chi         - component of node vector in eqw
%    psi         - component of node vector in eqw
%    meanlonM    - mean longitude                             rad
%    menalonNu   - true longitude                             rad
%    fr          - retrograde factor, neg if incl > 179 deg 1, -1
%
%  locals        :
%    none        -
%
%  coupling      :
%    none        -
%
%  references    :
%    vallado       2022, 110
%    chobotov            30
%
% [a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr] = rv2eq ( r, v );
% ----------------------------------------------------------------------------

function [a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr] = rv2eq ( r, v )

    constmath;
    constastro;

    % -------- convert to classical elements ----------------------
    [p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r, v);

    % -------- setup retrograde factor ----------------------------
    fr = 1.0;
    % ---------- set this so it is -1 only for orbits near 180 deg !! ---------
    if abs(incl - pi) < 0.0001
        fr = -1.0;
        incl*rad
    end

    % truelon = 0.0;
    % lonper = 0.0;
    % arglat = 0.0;

    coe = true;  % pick coe or vector approaches

    if (coe == true)
        if ( ecc < small )
            % ----------------  circular equatorial  ------------------
            if (incl<small) || ( abs(incl-pi)< small )
                argp = 0.0;
                omega= 0.0;
                nu   = truelon;
                m   = truelon;
            else
                % --------------  circular inclined  ------------------
                argp= 0.0;
                nu  = arglat;
                m  = arglat;
            end
        else
            % ---------------  elliptical equatorial  -----------------
            if ( ( incl<small) || (abs(incl-pi)<small) )
                argp = lonper;
                omega= 0.0;
            end
        end

        n  = sqrt(mu/(a*a*a));

        af = ecc * cos(fr*omega + argp);
        ag = ecc * sin(fr*omega + argp);

        if (fr > 0 )
            chi = tan(incl*0.5) * sin(omega);
            psi = tan(incl*0.5) * cos(omega);
        else
            chi = cot(incl*0.5) * sin(omega);
            psi = cot(incl*0.5) * cos(omega);
        end

        meanlonM = fr*omega + argp + m;
        meanlonM = rem(meanlonM+twopi, twopi);

        meanlonNu = fr*omega + argp + nu;
        meanlonNu = rem(meanlonNu+twopi, twopi);

        [eccanom, nu] = newtonm ( ecc, m );
        %         fprintf(1,'rv2eq F %11.7f  L %11.7f \n', (fr*omega + argp + eccanom)*rad, meanlonNu*rad);
        %         fprintf(1,'rv2eq F %11.7f  L %11.7f \n', (fr*omega + argp + eccanom)*rad, (fr*omega + argp + nu)*rad);

    else  % do with vectors
        magr = mag(r);
        magv = mag(v);

        a = 1.0 / (2.0/magr - magv^2/mu);
        n  = sqrt(mu/(a*a*a));

        wvec = cross(r,v) / mag(cross(r,v));
        chi = wvec(1)/(1.0 + fr*wvec(3));
        psi = -wvec(2)/(1.0 + fr*wvec(3));

        p0 = 1.0 / (1.0 + chi^2 + psi^2);
        fe = p0 * (1.0 - chi^2 + psi^2);  % 2nd one is minus???? no
        fq = p0 * 2.0 * chi * psi;
        fw = p0 * -2.0 * fr*chi;
        fvec = [fe, fq, fw];
        ge = p0 * 2.0 * fr*chi * psi;
        gq = p0 * fr*(1.0 + chi^2 - psi^2);
        gw = p0 * 2.0 * psi;
        gvec = [ge, gq, gw];
        we = wvec(1);
        wq = wvec(2);
        ww = wvec(3);
        we = p0 * 2.0 * chi;
        wq = p0 * -2.0 * psi;
        ww = p0 * fr*(1.0 - chi^2 - psi^2);
        wvec1 = [we, wq, ww];  % same

        % alt formulation - same as other, EXCEPT for retrograde (fr-1) so
        % do not use!!
        %         fe = 1.0 - we^2/(1.0 + ww);
        %         fq = -we*wq / (1.0 + ww);
        %         fw = -we;
        %         fvec = [fe, fq, fw];
        %         gvec = cross(wvec, fvec);

        evec = -r/magr + cross(v,cross(r,v)) / mu;

        ag = dot(evec, gvec);
        af = dot(evec, fvec);

        X = dot(r, fvec);
        Y = dot(r, gvec);

        temp = (af^2 + ag^2);
        if abs(temp) > 1.0
            temp = 1.0;
        end
        b = 1.0 / (1.0 + sqrt(1.0 - temp));

        sinF =  ag + ((1.0-ag^2*b)*Y - ag*af*b*X)/(a*sqrt(1.0 - temp));
        cosF =  af + ((1.0-af^2*b)*X - ag*af*b*Y)/(a*sqrt(1.0 - temp));
        F = atan2( sinF, cosF);

        %        fprintf(1,'rv2eq fe %11.7f %11.7f %11.7f ge %11.7f  %11.7f %11.7f \n X %11.7f Y %11.7f b %11.7f sF %11.7f cF %11.7f \n', fe, fq, fw, ge, gq, gw, X, Y, b, sinF, cosF);
        %        fprintf(1,'F = 316.20515  L = 13.61834  M = 288.88793 \n');

        sinZeta = ag / sqrt(temp);
        cosZeta = af / sqrt(temp);
        zeta = atan2( sinZeta, cosZeta );

        meanlonM = F + ag*cos(F) - af*sin(F);
        if meanlonM < 0.0
            meanlonM = twopi + meanlonM;
        end
        Eccanom = F - zeta;
        M = Eccanom - ecc*sin(Eccanom);
        if M < 0.0
            M = 2.0*pi + M;
        end
        %M = meanlonM - zeta;  % same

        sinL = ((1.0 - af^2*b)*sin(F) + ag*af*b*cos(F) - ag) / (1.0 - ag*sin(F) - af*cos(F));
        cosL = ((1.0 - ag^2*b)*cos(F) + ag*af*b*sin(F) - af) / (1.0 - ag*sin(F) - af*cos(F));
        L = atan2(sinL, cosL);

        meanlonNu = L;
        if meanlonNu < 0.0
            meanlonNu = twopi + meanlonNu;
        end
        nu = L - zeta;

        %       fprintf(1,'rv2eq F %11.7f  L %11.7f  %11.7f %11.7f %11.7f %11.7f  \n', F*rad, L*rad, X, Y, zeta*rad, M*rad);

    end

end